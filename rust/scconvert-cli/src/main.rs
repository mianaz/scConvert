//! scconvert — command-line front end.
//!
//! B1 scaffolding. Dispatches `convert`, `plan`, and `inspect` to
//! scconvert-core, which currently returns `NotYetImplemented` for
//! everything that does real IO. The binary exists so downstream packaging
//! (conda, Homebrew, `inst/bin/`) can wire up the *interface* before
//! format work lands in B2.

use scconvert_core::{ConvertOptions, Format};
use std::path::PathBuf;
use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("scconvert: {e:#}");
            // Map ConvertError codes into process exit codes >= 10 so
            // `0` stays reserved for success and `1` for CLI usage errors.
            if let Some(code) = e
                .downcast_ref::<scconvert_core::ConvertError>()
                .map(|c| c.code())
            {
                return ExitCode::from((code as u8).saturating_add(10));
            }
            ExitCode::from(1)
        }
    }
}

fn run(args: &[String]) -> anyhow::Result<()> {
    let Some((cmd, rest)) = args.split_first() else {
        print_usage();
        anyhow::bail!("missing subcommand");
    };

    match cmd.as_str() {
        "convert" => cmd_convert(rest),
        "plan" => cmd_plan(rest),
        "inspect" => cmd_inspect(rest),
        "--help" | "-h" | "help" => {
            print_usage();
            Ok(())
        }
        "--version" | "-V" => {
            println!("scconvert {}", env!("CARGO_PKG_VERSION"));
            Ok(())
        }
        other => {
            print_usage();
            anyhow::bail!("unknown subcommand: {other}");
        }
    }
}

fn cmd_convert(rest: &[String]) -> anyhow::Result<()> {
    let (src, dst, opts) = parse_src_dst_opts(rest)?;
    let report = scconvert_core::convert(&src, &dst, &opts)?;
    let json = report.to_json(true)?;
    println!("{json}");
    Ok(())
}

fn cmd_plan(rest: &[String]) -> anyhow::Result<()> {
    let (src, dst, opts) = parse_src_dst_opts(rest)?;
    let plan = scconvert_core::plan(&src, &dst, &opts)?;
    let json = serde_json::to_string_pretty(&plan)?;
    println!("{json}");
    Ok(())
}

fn cmd_inspect(rest: &[String]) -> anyhow::Result<()> {
    let [path] = rest else {
        anyhow::bail!("usage: scconvert inspect <path>");
    };
    let fmt = scconvert_core::detect_format(std::path::Path::new(path));
    #[derive(serde::Serialize)]
    struct Info<'a> {
        path: &'a str,
        format: Format,
    }
    let j = serde_json::to_string_pretty(&Info {
        path: path.as_str(),
        format: fmt,
    })?;
    println!("{j}");
    Ok(())
}

fn parse_src_dst_opts(rest: &[String]) -> anyhow::Result<(PathBuf, PathBuf, ConvertOptions)> {
    // Minimal flag parser for B1. Real arg parsing with `clap` can land
    // in B5 alongside `--loss-report`, `--max-memory`, `--threads`, etc.
    let mut positional: Vec<&String> = Vec::new();
    let mut opts = ConvertOptions::default();
    let mut i = 0;
    while i < rest.len() {
        match rest[i].as_str() {
            "--overwrite" => {
                opts.overwrite = true;
                i += 1;
            }
            "--threads" => {
                let v = rest
                    .get(i + 1)
                    .ok_or_else(|| anyhow::anyhow!("--threads needs a value"))?;
                opts.threads = v.parse()?;
                i += 2;
            }
            "--max-memory" => {
                let v = rest
                    .get(i + 1)
                    .ok_or_else(|| anyhow::anyhow!("--max-memory needs a value"))?;
                opts.max_memory_bytes = parse_bytes(v)?;
                i += 2;
            }
            "--compression" => {
                let v = rest
                    .get(i + 1)
                    .ok_or_else(|| anyhow::anyhow!("--compression needs a value"))?;
                opts.compression_level = Some(v.parse()?);
                i += 2;
            }
            other if other.starts_with("--") => {
                anyhow::bail!("unknown flag: {other}");
            }
            _ => {
                positional.push(&rest[i]);
                i += 1;
            }
        }
    }
    let [src, dst] = positional.as_slice() else {
        anyhow::bail!("usage: scconvert <convert|plan> <src> <dst> [flags]");
    };
    Ok((PathBuf::from(src), PathBuf::from(dst), opts))
}

fn parse_bytes(s: &str) -> anyhow::Result<u64> {
    // Accept plain bytes, `K`, `M`, `G` (binary) suffixes.
    let s = s.trim();
    let (num, mult) = if let Some(rest) = s.strip_suffix(['K', 'k']) {
        (rest, 1024u64)
    } else if let Some(rest) = s.strip_suffix(['M', 'm']) {
        (rest, 1024u64 * 1024)
    } else if let Some(rest) = s.strip_suffix(['G', 'g']) {
        (rest, 1024u64 * 1024 * 1024)
    } else {
        (s, 1)
    };
    let n: u64 = num.trim().parse()?;
    Ok(n * mult)
}

fn print_usage() {
    let exe = env!("CARGO_PKG_NAME");
    eprintln!(
        "\
{exe} {ver}

USAGE:
    scconvert convert <src> <dst> [--overwrite] [--threads N]
                                  [--max-memory 8G] [--compression 4]
    scconvert plan    <src> <dst>
    scconvert inspect <path>
    scconvert --version
    scconvert --help

B1 scaffolding. Real conversion paths land in B2 (h5ad<->h5seurat).
",
        ver = env!("CARGO_PKG_VERSION"),
    );
}
