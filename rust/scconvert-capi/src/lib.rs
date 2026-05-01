//! Stable C ABI over [`scconvert_core`].
//!
//! # Contract
//!
//! Every exported function returns an `i32` error code. `0` is success.
//! On error, callers receive a heap-allocated JSON string via an output
//! pointer, which they must free with [`sc_free_string`].
//!
//! Paths are passed as NUL-terminated UTF-8 `const char *`. The core
//! does not borrow them past the call boundary; callers may free their
//! buffers immediately afterwards.
//!
//! Options and results are JSON. This keeps the ABI small and stable
//! while letting the Rust side evolve structs freely behind serde.
//!
//! # Versioning
//!
//! Error codes are stable from 0.1 (see `ConvertError::code`). Adding a
//! new code is a minor bump; renumbering is a major bump.

use scconvert_core::{ConvertError, ConvertOptions};
use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int};
use std::path::Path;
use std::ptr;

/// ABI version. Bump the minor component on additive changes (new
/// functions, new error codes). Bump the major on any breaking change.
#[no_mangle]
pub extern "C" fn sc_abi_version() -> u32 {
    // 0x00_01_00_00 == 0.1.0
    0x0001_0000
}

/// Free a C string previously returned by this library.
///
/// # Safety
/// `s` must be a pointer previously returned by a `sc_*` function. Passing
/// `NULL` is a no-op. Double-free is UB.
#[no_mangle]
pub unsafe extern "C" fn sc_free_string(s: *mut c_char) {
    if !s.is_null() {
        // Retake ownership and drop.
        let _ = CString::from_raw(s);
    }
}

/// Convert `src` → `dst` using options encoded as a JSON string.
///
/// On success, writes a newly-allocated JSON fidelity report to
/// `out_report` and returns 0. On failure, writes a newly-allocated JSON
/// error blob to `out_error` and returns the [`ConvertError::code`]. The
/// caller owns the returned pointer and must free with
/// [`sc_free_string`]. Either pointer may be NULL to suppress that output.
///
/// # Safety
/// `src`, `dst`, and `options_json` must be NUL-terminated UTF-8.
/// `out_report` and `out_error`, if non-NULL, must point to valid
/// `*mut c_char` locations.
#[no_mangle]
pub unsafe extern "C" fn sc_convert(
    src: *const c_char,
    dst: *const c_char,
    options_json: *const c_char,
    out_report: *mut *mut c_char,
    out_error: *mut *mut c_char,
) -> c_int {
    clear_out(out_report);
    clear_out(out_error);

    let src = match read_path(src, "src") {
        Ok(p) => p,
        Err(e) => return emit_err(out_error, &e),
    };
    let dst = match read_path(dst, "dst") {
        Ok(p) => p,
        Err(e) => return emit_err(out_error, &e),
    };
    let opts = match read_options(options_json) {
        Ok(o) => o,
        Err(e) => return emit_err(out_error, &e),
    };

    match scconvert_core::convert(&src, &dst, &opts) {
        Ok(report) => match report.to_json(false) {
            Ok(j) => {
                emit_ok(out_report, &j);
                0
            }
            Err(e) => emit_err(out_error, &e),
        },
        Err(e) => emit_err(out_error, &e),
    }
}

/// Plan a conversion without executing. JSON-in, JSON-out.
///
/// # Safety
/// Same as [`sc_convert`].
#[no_mangle]
pub unsafe extern "C" fn sc_plan(
    src: *const c_char,
    dst: *const c_char,
    options_json: *const c_char,
    out_plan: *mut *mut c_char,
    out_error: *mut *mut c_char,
) -> c_int {
    clear_out(out_plan);
    clear_out(out_error);

    let src = match read_path(src, "src") {
        Ok(p) => p,
        Err(e) => return emit_err(out_error, &e),
    };
    let dst = match read_path(dst, "dst") {
        Ok(p) => p,
        Err(e) => return emit_err(out_error, &e),
    };
    let opts = match read_options(options_json) {
        Ok(o) => o,
        Err(e) => return emit_err(out_error, &e),
    };

    match scconvert_core::plan(&src, &dst, &opts) {
        Ok(plan) => match serde_json::to_string(&plan) {
            Ok(j) => {
                emit_ok(out_plan, &j);
                0
            }
            Err(e) => emit_err(out_error, &ConvertError::from(e)),
        },
        Err(e) => emit_err(out_error, &e),
    }
}

// ---------- helpers ----------

unsafe fn clear_out(out: *mut *mut c_char) {
    if !out.is_null() {
        *out = ptr::null_mut();
    }
}

fn read_path(raw: *const c_char, label: &str) -> Result<std::path::PathBuf, ConvertError> {
    if raw.is_null() {
        return Err(ConvertError::Other(format!("{} pointer is NULL", label)));
    }
    let cstr = unsafe { CStr::from_ptr(raw) };
    let s = cstr
        .to_str()
        .map_err(|_| ConvertError::Other(format!("{} is not valid UTF-8", label)))?;
    Ok(Path::new(s).to_path_buf())
}

fn read_options(raw: *const c_char) -> Result<ConvertOptions, ConvertError> {
    if raw.is_null() {
        return Ok(ConvertOptions::default());
    }
    let cstr = unsafe { CStr::from_ptr(raw) };
    let s = cstr
        .to_str()
        .map_err(|_| ConvertError::Other("options_json is not valid UTF-8".into()))?;
    if s.trim().is_empty() {
        return Ok(ConvertOptions::default());
    }
    Ok(serde_json::from_str(s)?)
}

unsafe fn emit_ok(out: *mut *mut c_char, s: &str) {
    if out.is_null() {
        return;
    }
    let Ok(cs) = CString::new(s) else { return };
    *out = cs.into_raw();
}

unsafe fn emit_err(out: *mut *mut c_char, e: &ConvertError) -> c_int {
    #[derive(serde::Serialize)]
    struct ErrBody<'a> {
        code: i32,
        message: &'a str,
    }
    let body = ErrBody {
        code: e.code(),
        message: &format!("{}", e),
    };
    let j = serde_json::to_string(&body)
        .unwrap_or_else(|_| r#"{"code":99,"message":"error serialisation failed"}"#.into());
    if !out.is_null() {
        if let Ok(cs) = CString::new(j) {
            *out = cs.into_raw();
        }
    }
    e.code()
}
