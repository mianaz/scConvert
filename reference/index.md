# Package index

## Package Information

Overview and utilities

- [`scConvert-package`](https://mianaz.github.io/scConvert/reference/scConvert-package.md)
  : scConvert: Universal Single-Cell Format Conversion for R

- [`scdisk-class`](https://mianaz.github.io/scConvert/reference/scdisk-class.md)
  [`scdisk`](https://mianaz.github.io/scConvert/reference/scdisk-class.md)
  : A disk-based object for single-cell analysis

- [`IsSCDisk()`](https://mianaz.github.io/scConvert/reference/IsSCDisk.md)
  : Does an R6 class inherit from scdisk

- [`GetSCDisk()`](https://mianaz.github.io/scConvert/reference/RegisterSCDisk.md)
  [`RegisterSCDisk()`](https://mianaz.github.io/scConvert/reference/RegisterSCDisk.md)
  :

  Get and Register `scdisk` Subclasses

## h5Seurat Format

Working with h5Seurat files - the native HDF5-based R/Seurat format

- [`writeH5Seurat()`](https://mianaz.github.io/scConvert/reference/writeH5Seurat.md)
  [`as.h5Seurat()`](https://mianaz.github.io/scConvert/reference/writeH5Seurat.md)
  : Save a Seurat object to an h5Seurat file
- [`readH5Seurat()`](https://mianaz.github.io/scConvert/reference/readH5Seurat.md)
  [`as.Seurat(`*`<h5Seurat>`*`)`](https://mianaz.github.io/scConvert/reference/readH5Seurat.md)
  : Load a Seurat object from an h5Seurat file
- [`h5Seurat-class`](https://mianaz.github.io/scConvert/reference/h5Seurat-class.md)
  [`h5Seurat`](https://mianaz.github.io/scConvert/reference/h5Seurat-class.md)
  : A class for connections to h5Seurat files
- [`Cells(`*`<h5Seurat>`*`)`](https://mianaz.github.io/scConvert/reference/h5Seurat-bindings.md)
  [`DefaultAssay(`*`<h5Seurat>`*`)`](https://mianaz.github.io/scConvert/reference/h5Seurat-bindings.md)
  [`` `DefaultAssay<-`( ``*`<h5Seurat>`*`)`](https://mianaz.github.io/scConvert/reference/h5Seurat-bindings.md)
  [`Idents(`*`<h5Seurat>`*`)`](https://mianaz.github.io/scConvert/reference/h5Seurat-bindings.md)
  [`IsGlobal(`*`<H5Group>`*`)`](https://mianaz.github.io/scConvert/reference/h5Seurat-bindings.md)
  [`Key(`*`<H5Group>`*`)`](https://mianaz.github.io/scConvert/reference/h5Seurat-bindings.md)
  [`Project(`*`<h5Seurat>`*`)`](https://mianaz.github.io/scConvert/reference/h5Seurat-bindings.md)
  [`` `Project<-`( ``*`<h5Seurat>`*`)`](https://mianaz.github.io/scConvert/reference/h5Seurat-bindings.md)
  [`Stdev(`*`<h5Seurat>`*`)`](https://mianaz.github.io/scConvert/reference/h5Seurat-bindings.md)
  : Seurat bindings for h5Seurat files
- [`scConnect()`](https://mianaz.github.io/scConvert/reference/scConnect.md)
  : scConnect to a single-cell HDF5 dataset
- [`DefaultAssay(`*`<h5SI>`*`)`](https://mianaz.github.io/scConvert/reference/h5SI.md)
  [`print(`*`<h5SI>`*`)`](https://mianaz.github.io/scConvert/reference/h5SI.md)
  : Tools for handling h5Seurat indexes

## AnnData (h5ad) Format

Converting between Seurat and Python AnnData format

- [`H5ADToH5Seurat()`](https://mianaz.github.io/scConvert/reference/H5ADToH5Seurat.md)
  : Convert AnnData/H5AD files to h5Seurat files
- [`H5SeuratToH5AD()`](https://mianaz.github.io/scConvert/reference/H5SeuratToH5AD.md)
  : Convert h5Seurat files to H5AD files
- [`writeH5AD()`](https://mianaz.github.io/scConvert/reference/writeH5AD.md)
  : Direct Seurat to H5AD Conversion
- [`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md)
  : Load an AnnData H5AD file as a Seurat object
- [`readH5AD_obs()`](https://mianaz.github.io/scConvert/reference/readH5AD_obs.md)
  : Read observation metadata from an H5AD file
- [`readH5AD_obsm()`](https://mianaz.github.io/scConvert/reference/readH5AD_obsm.md)
  : Read observation embeddings from an H5AD file
- [`scLoadMeta()`](https://mianaz.github.io/scConvert/reference/scLoadMeta.md)
  : Load deferred h5ad components into a Seurat object

## MuData (h5mu) Format

Working with MuData multimodal HDF5 files

- [`readH5MU()`](https://mianaz.github.io/scConvert/reference/readH5MU.md)
  : Load a MuData H5MU file as a Seurat object
- [`writeH5MU()`](https://mianaz.github.io/scConvert/reference/writeH5MU.md)
  : Write a Seurat object to h5mu format
- [`as.h5mu()`](https://mianaz.github.io/scConvert/reference/as.h5mu.md)
  : Save a Seurat object to h5mu format
- [`H5MUToH5Seurat()`](https://mianaz.github.io/scConvert/reference/H5MUToH5Seurat.md)
  : Convert an h5mu file to h5seurat format
- [`H5SeuratToH5MU()`](https://mianaz.github.io/scConvert/reference/H5SeuratToH5MU.md)
  : Convert an h5seurat file to h5mu format
- [`H5MUToH5AD()`](https://mianaz.github.io/scConvert/reference/H5MUToH5AD.md)
  : Convert an h5mu file to h5ad format
- [`H5ADToH5MU()`](https://mianaz.github.io/scConvert/reference/H5ADToH5MU.md)
  : Convert an h5ad file to h5mu format
- [`H5MUToLoom()`](https://mianaz.github.io/scConvert/reference/H5MUToLoom.md)
  : Convert an h5mu file to loom format
- [`LoomToH5MU()`](https://mianaz.github.io/scConvert/reference/LoomToH5MU.md)
  : Convert a loom file to h5mu format
- [`H5MUToZarr()`](https://mianaz.github.io/scConvert/reference/H5MUToZarr.md)
  : Convert an h5mu file to zarr format
- [`ZarrToH5MU()`](https://mianaz.github.io/scConvert/reference/ZarrToH5MU.md)
  : Convert a zarr store to h5mu format
- [`H5MUToSOMA()`](https://mianaz.github.io/scConvert/reference/H5MUToSOMA.md)
  : Convert H5MU (MuData) to SOMA format
- [`SOMAToH5MU()`](https://mianaz.github.io/scConvert/reference/SOMAToH5MU.md)
  : Convert SOMA to H5MU (MuData) format

## Loom Format

Loading and saving Loom files for interoperability with loompy

- [`readLoom()`](https://mianaz.github.io/scConvert/reference/readLoom.md)
  [`as.Seurat(`*`<loom>`*`)`](https://mianaz.github.io/scConvert/reference/readLoom.md)
  : Load a Loom file as a Seurat object
- [`writeLoom()`](https://mianaz.github.io/scConvert/reference/writeLoom.md)
  [`as.loom()`](https://mianaz.github.io/scConvert/reference/writeLoom.md)
  : Save a Seurat object to a Loom file
- [`loom-class`](https://mianaz.github.io/scConvert/reference/loom-class.md)
  [`loom`](https://mianaz.github.io/scConvert/reference/loom-class.md) :
  A class for connections to loom files
- [`DefaultAssay(`*`<loom>`*`)`](https://mianaz.github.io/scConvert/reference/loom-bindings.md)
  [`dim(`*`<loom>`*`)`](https://mianaz.github.io/scConvert/reference/loom-bindings.md)
  : Seurat binding for loom files
- [`LoadLoom0.1()`](https://mianaz.github.io/scConvert/reference/LoomLoading.md)
  [`LoadLoom3.0()`](https://mianaz.github.io/scConvert/reference/LoomLoading.md)
  : Loom-file Loading
- [`LoomValidate0.1()`](https://mianaz.github.io/scConvert/reference/ValidateLoom.md)
  [`LoomValidate3.0.0()`](https://mianaz.github.io/scConvert/reference/ValidateLoom.md)
  : Validate Loom Files
- [`LoomToH5Seurat()`](https://mianaz.github.io/scConvert/reference/LoomToH5Seurat.md)
  : Convert a loom file to h5seurat format
- [`H5SeuratToLoom()`](https://mianaz.github.io/scConvert/reference/H5SeuratToLoom.md)
  : Convert an h5seurat file to loom format
- [`LoomToH5AD()`](https://mianaz.github.io/scConvert/reference/LoomToH5AD.md)
  : Convert a loom file to h5ad format
- [`H5ADToLoom()`](https://mianaz.github.io/scConvert/reference/H5ADToLoom.md)
  : Convert an h5ad file to loom format
- [`LoomToZarr()`](https://mianaz.github.io/scConvert/reference/LoomToZarr.md)
  : Convert a loom file to zarr format
- [`ZarrToLoom()`](https://mianaz.github.io/scConvert/reference/ZarrToLoom.md)
  : Convert a zarr store to loom format
- [`LoomToSOMA()`](https://mianaz.github.io/scConvert/reference/LoomToSOMA.md)
  : Convert Loom to SOMA format
- [`SOMAToLoom()`](https://mianaz.github.io/scConvert/reference/SOMAToLoom.md)
  : Convert SOMA to Loom format

## Zarr Format

Loading, saving, and converting Zarr-based AnnData stores

- [`readZarr()`](https://mianaz.github.io/scConvert/reference/readZarr.md)
  : Load an AnnData Zarr store as a Seurat object
- [`writeZarr()`](https://mianaz.github.io/scConvert/reference/writeZarr.md)
  : Save a Seurat object as an AnnData Zarr store
- [`H5ADToZarr()`](https://mianaz.github.io/scConvert/reference/H5ADToZarr.md)
  : Convert an h5ad file to zarr format
- [`ZarrToH5AD()`](https://mianaz.github.io/scConvert/reference/ZarrToH5AD.md)
  : Convert a zarr store to h5ad format
- [`H5SeuratToZarr()`](https://mianaz.github.io/scConvert/reference/H5SeuratToZarr.md)
  : Convert an h5seurat file to zarr (AnnData) format
- [`ZarrToH5Seurat()`](https://mianaz.github.io/scConvert/reference/ZarrToH5Seurat.md)
  : Convert a zarr (AnnData) store to h5seurat format
- [`ZarrToSOMA()`](https://mianaz.github.io/scConvert/reference/ZarrToSOMA.md)
  : Convert Zarr to SOMA format
- [`SOMAToZarr()`](https://mianaz.github.io/scConvert/reference/SOMAToZarr.md)
  : Convert SOMA to Zarr format

## TileDB-SOMA Format

Reading and writing TileDB-SOMA for CELLxGENE Census interoperability

- [`readSOMA()`](https://mianaz.github.io/scConvert/reference/readSOMA.md)
  : Read a TileDB-SOMA experiment as a Seurat object
- [`writeSOMA()`](https://mianaz.github.io/scConvert/reference/writeSOMA.md)
  : Write a Seurat object to TileDB-SOMA format
- [`H5ADToSOMA()`](https://mianaz.github.io/scConvert/reference/H5ADToSOMA.md)
  : Convert h5ad to SOMA format
- [`SOMAToH5AD()`](https://mianaz.github.io/scConvert/reference/SOMAToH5AD.md)
  : Convert SOMA to h5ad format
- [`H5SeuratToSOMA()`](https://mianaz.github.io/scConvert/reference/H5SeuratToSOMA.md)
  : Convert h5Seurat to SOMA format
- [`SOMAToH5Seurat()`](https://mianaz.github.io/scConvert/reference/SOMAToH5Seurat.md)
  : Convert SOMA to h5Seurat format

## Spatial Data Conversion

Converting spatial data between Seurat, AnnData, and SpatialData formats

- [`H5ADSpatialToSeurat()`](https://mianaz.github.io/scConvert/reference/H5ADSpatialToSeurat.md)
  : Convert spatial coordinates from h5ad to Seurat format
- [`SeuratSpatialToH5AD()`](https://mianaz.github.io/scConvert/reference/SeuratSpatialToH5AD.md)
  : Convert Seurat spatial data to h5ad format
- [`readSpatialData()`](https://mianaz.github.io/scConvert/reference/readSpatialData.md)
  : Read a SpatialData zarr store into a Seurat object
- [`writeSpatialData()`](https://mianaz.github.io/scConvert/reference/writeSpatialData.md)
  : Write a Seurat object to SpatialData zarr format
- [`SpatialDataToH5AD()`](https://mianaz.github.io/scConvert/reference/SpatialDataToH5AD.md)
  : Convert a SpatialData zarr store to h5ad format
- [`H5ADToSpatialData()`](https://mianaz.github.io/scConvert/reference/H5ADToSpatialData.md)
  : Convert an h5ad file to SpatialData zarr format
- [`SpatialDataToH5Seurat()`](https://mianaz.github.io/scConvert/reference/SpatialDataToH5Seurat.md)
  : Convert a SpatialData zarr store to h5Seurat format
- [`H5SeuratToSpatialData()`](https://mianaz.github.io/scConvert/reference/H5SeuratToSpatialData.md)
  : Convert an h5Seurat file to SpatialData zarr format
- [`SpatialDataToZarr()`](https://mianaz.github.io/scConvert/reference/SpatialDataToZarr.md)
  : Convert a SpatialData zarr store to a standard anndata zarr store
- [`ZarrToSpatialData()`](https://mianaz.github.io/scConvert/reference/ZarrToSpatialData.md)
  : Wrap a standard anndata zarr store as a SpatialData store
- [`WriteFOVToH5AD()`](https://mianaz.github.io/scConvert/reference/WriteFOVToH5AD.md)
  : Write FOV segmentation and molecules into an open h5ad library group
- [`ReadFOVFromH5AD()`](https://mianaz.github.io/scConvert/reference/ReadFOVFromH5AD.md)
  : Read FOV structure back from an h5ad library group

## Vendor Spatial Formats

Native readers for upstream vendor formats (no Python required)

- [`LoadStereoSeqGef()`](https://mianaz.github.io/scConvert/reference/LoadStereoSeqGef.md)
  : Load a Stereo-seq GEF file into a Seurat object
- [`LoadCosMx()`](https://mianaz.github.io/scConvert/reference/LoadCosMx.md)
  : Load a NanoString CosMx SMI bundle into a Seurat object
- [`LoadXenium()`](https://mianaz.github.io/scConvert/reference/LoadXenium.md)
  : Load a 10x Xenium output bundle into a Seurat object

## Format Conversion

General conversion utilities between all supported formats

- [`scConvert()`](https://mianaz.github.io/scConvert/reference/scConvert.md)
  : Convert single-cell datasets between formats
- [`FileType()`](https://mianaz.github.io/scConvert/reference/FileType.md)
  : Determine a filetype based on its extension
- [`scConvert_cli()`](https://mianaz.github.io/scConvert/reference/scConvert_cli.md)
  : Run the scconvert CLI for file-to-file conversion

## Reading & Writing

Low-level HDF5 reading and writing functions

- [`as.array(`*`<H5D>`*`)`](https://mianaz.github.io/scConvert/reference/ReadH5.md)
  [`as.data.frame(`*`<H5D>`*`)`](https://mianaz.github.io/scConvert/reference/ReadH5.md)
  [`as.data.frame(`*`<H5Group>`*`)`](https://mianaz.github.io/scConvert/reference/ReadH5.md)
  [`as.list(`*`<H5D>`*`)`](https://mianaz.github.io/scConvert/reference/ReadH5.md)
  [`as.list(`*`<H5Group>`*`)`](https://mianaz.github.io/scConvert/reference/ReadH5.md)
  [`as.matrix(`*`<H5D>`*`)`](https://mianaz.github.io/scConvert/reference/ReadH5.md)
  [`as.matrix(`*`<H5Group>`*`)`](https://mianaz.github.io/scConvert/reference/ReadH5.md)
  : Read HDF5 Files
- [`as.sparse(`*`<H5D>`*`)`](https://mianaz.github.io/scConvert/reference/as.sparse.md)
  [`as.sparse(`*`<H5Group>`*`)`](https://mianaz.github.io/scConvert/reference/as.sparse.md)
  : Convert an HDF5 dataset to a sparse matrix
- [`WriteH5Group()`](https://mianaz.github.io/scConvert/reference/WriteH5Group.md)
  : Write data to an HDF5 group
- [`BasicWrite()`](https://mianaz.github.io/scConvert/reference/BasicWrite.md)
  : Write lists and other data to an HDF5 dataset
- [`ImageWrite()`](https://mianaz.github.io/scConvert/reference/ImageWrite.md)
  : Write a SpatialImage object to an HDF5 dataset
- [`SparseWrite()`](https://mianaz.github.io/scConvert/reference/SparseWrite.md)
  : Write a sparse matrix to an HDF5 dataset
- [`AttrExists()`](https://mianaz.github.io/scConvert/reference/H5Exists.md)
  [`Exists()`](https://mianaz.github.io/scConvert/reference/H5Exists.md)
  : Check to see if a dataset, group, or attribute exists in an HDF5
  file, group, or dataset
- [`H5Path()`](https://mianaz.github.io/scConvert/reference/H5Path.md) :
  Create an HDF5 object path

## Object Operations

Utilities for working with HDF5 objects and attributes

- [`GetAssays()`](https://mianaz.github.io/scConvert/reference/GetObject.md)
  [`GetCommands()`](https://mianaz.github.io/scConvert/reference/GetObject.md)
  [`GetDimReducs()`](https://mianaz.github.io/scConvert/reference/GetObject.md)
  [`GetGraphs()`](https://mianaz.github.io/scConvert/reference/GetObject.md)
  [`GetImages()`](https://mianaz.github.io/scConvert/reference/GetObject.md)
  [`GetNeighbors()`](https://mianaz.github.io/scConvert/reference/GetObject.md)
  : Figure out which objects to load from an h5Seurat file
- [`AssembleAssay()`](https://mianaz.github.io/scConvert/reference/AssembleObject.md)
  [`AssembleDimReduc()`](https://mianaz.github.io/scConvert/reference/AssembleObject.md)
  [`AssembleGraph()`](https://mianaz.github.io/scConvert/reference/AssembleObject.md)
  [`AssembleImage()`](https://mianaz.github.io/scConvert/reference/AssembleObject.md)
  [`AssembleNeighbor()`](https://mianaz.github.io/scConvert/reference/AssembleObject.md)
  [`AssembleSeuratCommand()`](https://mianaz.github.io/scConvert/reference/AssembleObject.md)
  : Assemble an object from an h5Seurat file
- [`scTranspose()`](https://mianaz.github.io/scConvert/reference/scTranspose.md)
  : scTranspose a matrix
- [`PadMatrix()`](https://mianaz.github.io/scConvert/reference/PadMatrix.md)
  : Pad a matrix

## Testing & Validation

Functions for testing and validating HDF5 files and Seurat objects

- [`IsDataFrame(`*`<H5D>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  [`IsDataFrame(`*`<H5Group>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  [`IsFactor(`*`<H5D>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  [`IsFactor(`*`<H5Group>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  [`IsList(`*`<H5Group>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  [`IsLogical(`*`<H5D>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  [`IsMatrix(`*`<H5D>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  [`IsMatrix(`*`<H5Group>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  : Test HDF5 datasets and groups to see what kind of data they are
- [`IsDataFrame()`](https://mianaz.github.io/scConvert/reference/TestObject.md)
  [`IsFactor()`](https://mianaz.github.io/scConvert/reference/TestObject.md)
  [`IsList()`](https://mianaz.github.io/scConvert/reference/TestObject.md)
  [`IsLogical()`](https://mianaz.github.io/scConvert/reference/TestObject.md)
  [`IsMatrix()`](https://mianaz.github.io/scConvert/reference/TestObject.md)
  : Test an object's class
- [`CheckMatrix()`](https://mianaz.github.io/scConvert/reference/CheckMatrix.md)
  : Check that a dataset is a proper loom matrix
- [`IsMatrixEmpty()`](https://mianaz.github.io/scConvert/reference/IsMatrixEmpty.md)
  : Check to see if a matrix is empty

## Data Type Utilities

Utilities for handling data types and format conversions

- [`GuessDType()`](https://mianaz.github.io/scConvert/reference/GuessDType.md)
  : Guess an HDF5 Datatype
- [`IsDType()`](https://mianaz.github.io/scConvert/reference/IsDType.md)
  : Check the datatype of an HDF5 dataset
- [`StringType()`](https://mianaz.github.io/scConvert/reference/StringType.md)
  : Generate an HDF5 string dtype
- [`BoolToInt()`](https://mianaz.github.io/scConvert/reference/BoolToInt.md)
  : Convert a logical to an integer
- [`Dims()`](https://mianaz.github.io/scConvert/reference/Dims.md) : Get
  the dimensions of an HDF5 dataset or sparse matrix
- [`GetClass()`](https://mianaz.github.io/scConvert/reference/GetClass.md)
  : Get a class string with package information
- [`GetMargin()`](https://mianaz.github.io/scConvert/reference/GetMargin.md)
  : Determine the margin to use for a dataset
- [`GetParent()`](https://mianaz.github.io/scConvert/reference/GetParent.md)
  : Get the parent of an HDF5 dataset or group
- [`IndexToPointer()`](https://mianaz.github.io/scConvert/reference/SparsePointers.md)
  [`PointerToIndex()`](https://mianaz.github.io/scConvert/reference/SparsePointers.md)
  : Convert sparse matrix pointers to indices and vice versa
- [`ClosestVersion()`](https://mianaz.github.io/scConvert/reference/ClosestVersion.md)
  : Find the closest version
- [`ChunkPoints()`](https://mianaz.github.io/scConvert/reference/ChunkPoints.md)
  : Generate chunk points
- [`MakeSpace()`](https://mianaz.github.io/scConvert/reference/MakeSpace.md)
  : Make a space

## Metadata & Attributes

Working with HDF5 attributes and metadata

- [`PadNames()`](https://mianaz.github.io/scConvert/reference/PadNames.md)
  : Add names for unnamed or partially named objects
- [`RandomName()`](https://mianaz.github.io/scConvert/reference/RandomName.md)
  : Generate a random string of characters
- [`FormatTime()`](https://mianaz.github.io/scConvert/reference/Timestamp.md)
  [`Timestamp()`](https://mianaz.github.io/scConvert/reference/Timestamp.md)
  [`TSFormats()`](https://mianaz.github.io/scConvert/reference/Timestamp.md)
  : Create and work with timestamps
- [`UpdateKey()`](https://mianaz.github.io/scConvert/reference/UpdateKey.md)
  : Update a Seurat key
- [`UpdateSlots()`](https://mianaz.github.io/scConvert/reference/UpdateSlots.md)
  : Update slots in an object
- [`FixFeatures()`](https://mianaz.github.io/scConvert/reference/FixFeatures.md)
  : Fix Feature Names
- [`CompoundToGroup()`](https://mianaz.github.io/scConvert/reference/CompoundToGroup.md)
  : Convert compound dataset to group
- [`Writeable()`](https://mianaz.github.io/scConvert/reference/Writeable.md)
  : Is an HDF5 file or group writeable
- [`WriteMode()`](https://mianaz.github.io/scConvert/reference/WriteMode.md)
  : Get the proper HDF5 connection mode for writing depending on
  overwrite status

## Internal

Internal functions and methods

- [`AddAnndataEncoding()`](https://mianaz.github.io/scConvert/reference/AddAnndataEncoding.md)
  : Add AnnData encoding attributes to an HDF5 dataset or group

- [`AddGenericSpatialData()`](https://mianaz.github.io/scConvert/reference/AddGenericSpatialData.md)
  : Add generic spatial data to Seurat object

- [`AddSlideSeqSpatialData()`](https://mianaz.github.io/scConvert/reference/AddSlideSeqSpatialData.md)
  : Add SlideSeq spatial data to Seurat object

- [`AddVisiumSpatialData()`](https://mianaz.github.io/scConvert/reference/AddVisiumSpatialData.md)
  : Add Visium spatial data to Seurat object

- [`AnnDataLayerToSeurat()`](https://mianaz.github.io/scConvert/reference/AnnDataLayerToSeurat.md)
  : Map AnnData layer name to Seurat layer name

- [`AnnDataReductionKey()`](https://mianaz.github.io/scConvert/reference/AnnDataReductionKey.md)
  : Map AnnData reduction name to Seurat key

- [`AssembleAssay()`](https://mianaz.github.io/scConvert/reference/AssembleObject.md)
  [`AssembleDimReduc()`](https://mianaz.github.io/scConvert/reference/AssembleObject.md)
  [`AssembleGraph()`](https://mianaz.github.io/scConvert/reference/AssembleObject.md)
  [`AssembleImage()`](https://mianaz.github.io/scConvert/reference/AssembleObject.md)
  [`AssembleNeighbor()`](https://mianaz.github.io/scConvert/reference/AssembleObject.md)
  [`AssembleSeuratCommand()`](https://mianaz.github.io/scConvert/reference/AssembleObject.md)
  : Assemble an object from an h5Seurat file

- [`BasicWrite()`](https://mianaz.github.io/scConvert/reference/BasicWrite.md)
  : Write lists and other data to an HDF5 dataset

- [`BoolToInt()`](https://mianaz.github.io/scConvert/reference/BoolToInt.md)
  : Convert a logical to an integer

- [`CachedGuessDType()`](https://mianaz.github.io/scConvert/reference/CachedGuessDType.md)
  : Cached version of GuessDType for known string constants

- [`CachedUtf8Type()`](https://mianaz.github.io/scConvert/reference/CachedUtf8Type.md)
  : Get cached UTF-8 string dtype

- [`CheckMatrix()`](https://mianaz.github.io/scConvert/reference/CheckMatrix.md)
  : Check that a dataset is a proper loom matrix

- [`ChunkPoints()`](https://mianaz.github.io/scConvert/reference/ChunkPoints.md)
  : Generate chunk points

- [`ClosestVersion()`](https://mianaz.github.io/scConvert/reference/ClosestVersion.md)
  : Find the closest version

- [`CompoundToGroup()`](https://mianaz.github.io/scConvert/reference/CompoundToGroup.md)
  : Convert compound dataset to group

- [`ConvertBPCellsMatrix()`](https://mianaz.github.io/scConvert/reference/ConvertBPCellsMatrix.md)
  : Convert BPCells matrix objects to dgCMatrix

- [`ConvertFOVToH5AD()`](https://mianaz.github.io/scConvert/reference/ConvertFOVToH5AD.md)
  : Convert FOV object to an h5ad-compatible list structure

- [`ConvertH5MUSpatialToSeurat()`](https://mianaz.github.io/scConvert/reference/ConvertH5MUSpatialToSeurat.md)
  : Convert h5mu spatial data to Seurat format for a specific modality

- [`CopyH5MatrixData()`](https://mianaz.github.io/scConvert/reference/CopyH5MatrixData.md)
  : Copy HDF5 matrix data handling nested paths

- [`CreateCachedExistsChecker()`](https://mianaz.github.io/scConvert/reference/CreateCachedExistsChecker.md)
  : Create a cached HDF5 existence checker

- [`DecodeCategorical()`](https://mianaz.github.io/scConvert/reference/DecodeCategorical.md)
  : Decode AnnData categorical encoding to R factor

- [`DeconstructSparseCSR()`](https://mianaz.github.io/scConvert/reference/DeconstructSparseCSR.md)
  : Deconstruct dgCMatrix to AnnData CSR components

- [`DetectMultiLibrary()`](https://mianaz.github.io/scConvert/reference/DetectMultiLibrary.md)
  : Detect if h5Seurat file contains multi-library spatial data

- [`DetectSpatialTechnology()`](https://mianaz.github.io/scConvert/reference/DetectSpatialTechnology.md)
  : Detect spatial technology from data structure

- [`Dims()`](https://mianaz.github.io/scConvert/reference/Dims.md) : Get
  the dimensions of an HDF5 dataset or sparse matrix

- [`DirectSeuratToH5AD()`](https://mianaz.github.io/scConvert/reference/DirectSeuratToH5AD.md)
  : Direct Seurat-to-h5ad write (single-pass, no intermediate h5Seurat)

- [`EncodeCategorical()`](https://mianaz.github.io/scConvert/reference/EncodeCategorical.md)
  : Encode R factor as AnnData categorical

- [`ExtractFOVCentroids()`](https://mianaz.github.io/scConvert/reference/ExtractFOVCentroids.md)
  : Extract centroid coordinates from an FOV object

- [`ExtractFOVMolecules()`](https://mianaz.github.io/scConvert/reference/ExtractFOVMolecules.md)
  : Extract molecule coordinates from an FOV object

- [`ExtractFOVSegmentation()`](https://mianaz.github.io/scConvert/reference/ExtractFOVSegmentation.md)
  : Extract segmentation boundaries from an FOV object

- [`FileType()`](https://mianaz.github.io/scConvert/reference/FileType.md)
  : Determine a filetype based on its extension

- [`FixFeatures()`](https://mianaz.github.io/scConvert/reference/FixFeatures.md)
  : Fix Feature Names

- [`GetAssayDataCompat()`](https://mianaz.github.io/scConvert/reference/GetAssayDataCompat.md)
  : Get assay data with V4/V5 compatibility

- [`GetAssayLayersCompat()`](https://mianaz.github.io/scConvert/reference/GetAssayLayersCompat.md)
  : Get available layers or slots from assay

- [`GetClass()`](https://mianaz.github.io/scConvert/reference/GetClass.md)
  : Get a class string with package information

- [`GetCompressionLevel()`](https://mianaz.github.io/scConvert/reference/GetCompressionLevel.md)
  : Get compression level for HDF5 dataset creation

- [`GetDefaultAssayToModalityMapping()`](https://mianaz.github.io/scConvert/reference/GetDefaultAssayToModalityMapping.md)
  : Get default assay to modality name mapping

- [`GetDefaultModalityMapping()`](https://mianaz.github.io/scConvert/reference/GetDefaultModalityMapping.md)
  : Get default modality to assay name mapping

- [`GetFeaturesV5Safe()`](https://mianaz.github.io/scConvert/reference/GetFeaturesV5Safe.md)
  : Get feature names with V5 fallback

- [`GetLayerPath()`](https://mianaz.github.io/scConvert/reference/GetLayerPath.md)
  : Get Layer Path

- [`GetMargin()`](https://mianaz.github.io/scConvert/reference/GetMargin.md)
  : Determine the margin to use for a dataset

- [`GetAssays()`](https://mianaz.github.io/scConvert/reference/GetObject.md)
  [`GetCommands()`](https://mianaz.github.io/scConvert/reference/GetObject.md)
  [`GetDimReducs()`](https://mianaz.github.io/scConvert/reference/GetObject.md)
  [`GetGraphs()`](https://mianaz.github.io/scConvert/reference/GetObject.md)
  [`GetImages()`](https://mianaz.github.io/scConvert/reference/GetObject.md)
  [`GetNeighbors()`](https://mianaz.github.io/scConvert/reference/GetObject.md)
  : Figure out which objects to load from an h5Seurat file

- [`GetParent()`](https://mianaz.github.io/scConvert/reference/GetParent.md)
  : Get the parent of an HDF5 dataset or group

- [`GetScaleFactors()`](https://mianaz.github.io/scConvert/reference/GetScaleFactors.md)
  : Get scale factors from Visium object

- [`GetStringType()`](https://mianaz.github.io/scConvert/reference/GetStringType.md)
  : Check if an HDF5 dataset is a factor

- [`GuessDType()`](https://mianaz.github.io/scConvert/reference/GuessDType.md)
  : Guess an HDF5 Datatype

- [`H5ADToH5Seurat()`](https://mianaz.github.io/scConvert/reference/H5ADToH5Seurat.md)
  : Convert AnnData/H5AD files to h5Seurat files

- [`AttrExists()`](https://mianaz.github.io/scConvert/reference/H5Exists.md)
  [`Exists()`](https://mianaz.github.io/scConvert/reference/H5Exists.md)
  : Check to see if a dataset, group, or attribute exists in an HDF5
  file, group, or dataset

- [`H5Path()`](https://mianaz.github.io/scConvert/reference/H5Path.md) :
  Create an HDF5 object path

- [`H5SeuratToH5AD()`](https://mianaz.github.io/scConvert/reference/H5SeuratToH5AD.md)
  : Convert h5Seurat files to H5AD files

- [`ImageWrite()`](https://mianaz.github.io/scConvert/reference/ImageWrite.md)
  : Write a SpatialImage object to an HDF5 dataset

- [`IsAssay5()`](https://mianaz.github.io/scConvert/reference/IsAssay5.md)
  : Check if object is Assay5 (V5 format)

- [`IsCompound()`](https://mianaz.github.io/scConvert/reference/IsCompound.md)
  : Check for compound type in HDF5 dataset

- [`IsDType()`](https://mianaz.github.io/scConvert/reference/IsDType.md)
  : Check the datatype of an HDF5 dataset

- [`IsEffectively1D()`](https://mianaz.github.io/scConvert/reference/IsEffectively1D.md)
  : Check if H5D dataset is effectively 1D

- [`IsMatrixEmpty()`](https://mianaz.github.io/scConvert/reference/IsMatrixEmpty.md)
  : Check to see if a matrix is empty

- [`IsSCDisk()`](https://mianaz.github.io/scConvert/reference/IsSCDisk.md)
  : Does an R6 class inherit from scdisk

- [`IsV5Structure()`](https://mianaz.github.io/scConvert/reference/IsV5Structure.md)
  : Check if h5Seurat file uses V5 structure

- [`ListFormats()`](https://mianaz.github.io/scConvert/reference/ListFormats.md)
  : List all registered format conversions

- [`LoadLoom0.1()`](https://mianaz.github.io/scConvert/reference/LoomLoading.md)
  [`LoadLoom3.0()`](https://mianaz.github.io/scConvert/reference/LoomLoading.md)
  : Loom-file Loading

- [`MakeSpace()`](https://mianaz.github.io/scConvert/reference/MakeSpace.md)
  : Make a space

- [`MigrateV4ToV5()`](https://mianaz.github.io/scConvert/reference/MigrateV4ToV5.md)
  : Migrate V4 to V5 Structure

- [`PB()`](https://mianaz.github.io/scConvert/reference/PB.md) : Create
  a progress bar

- [`PadMatrix()`](https://mianaz.github.io/scConvert/reference/PadMatrix.md)
  : Pad a matrix

- [`PadNames()`](https://mianaz.github.io/scConvert/reference/PadNames.md)
  : Add names for unnamed or partially named objects

- [`PrepareSpatialForH5MU()`](https://mianaz.github.io/scConvert/reference/PrepareSpatialForH5MU.md)
  : Prepare spatial data for h5mu export

- [`RandomName()`](https://mianaz.github.io/scConvert/reference/RandomName.md)
  : Generate a random string of characters

- [`ReadH5ImageDataset()`](https://mianaz.github.io/scConvert/reference/ReadH5ImageDataset.md)
  : Read HDF5 image dataset with proper dimension handling

- [`ReadLibraryIds()`](https://mianaz.github.io/scConvert/reference/ReadLibraryIds.md)
  : Read library IDs from h5Seurat metadata

- [`ReadSparseMatrix()`](https://mianaz.github.io/scConvert/reference/ReadSparseMatrix.md)
  : Read sparse matrix from HDF5

- [`ReadV5Layer()`](https://mianaz.github.io/scConvert/reference/ReadV5Layer.md)
  : Read V5 Layer Data

- [`ReconstructSparseCSR()`](https://mianaz.github.io/scConvert/reference/ReconstructSparseCSR.md)
  : Reconstruct sparse CSR/CSC matrix from array components

- [`RegisterDirectPath()`](https://mianaz.github.io/scConvert/reference/RegisterDirectPath.md)
  : Register a direct conversion path between two file formats

- [`RegisterFormat()`](https://mianaz.github.io/scConvert/reference/RegisterFormat.md)
  : Register a file format with its load and save functions

- [`RequiresS4Reconstruction()`](https://mianaz.github.io/scConvert/reference/RequiresS4Reconstruction.md)
  : Check if assay requires S4 slot reconstruction

- [`ResolveNestedH5Path()`](https://mianaz.github.io/scConvert/reference/ResolveNestedH5Path.md)
  : Resolve a nested HDF5 path to its object

- [`RestoreSpatialFromH5MU()`](https://mianaz.github.io/scConvert/reference/RestoreSpatialFromH5MU.md)
  : Restore spatial data from h5mu file to Seurat object

- [`SafeH5DRead()`](https://mianaz.github.io/scConvert/reference/SafeH5DRead.md)
  : Safely read an H5D dataset

- [`SafeH5GroupToList()`](https://mianaz.github.io/scConvert/reference/SafeH5GroupToList.md)
  : Safely convert H5Group to list, handling 3D+ arrays

- [`SafeSetLayerData()`](https://mianaz.github.io/scConvert/reference/SafeSetLayerData.md)
  : Safely set layer data in a V5 Assay

- [`SanitizeImageKey()`](https://mianaz.github.io/scConvert/reference/SanitizeImageKey.md)
  : Sanitize string to valid Seurat image key

- [`ScalarSpace()`](https://mianaz.github.io/scConvert/reference/ScalarSpace.md)
  : Get cached scalar HDF5 dataspace

- [`SetAssayDataCompat()`](https://mianaz.github.io/scConvert/reference/SetAssayDataCompat.md)
  : Set assay data with V4/V5 compatibility

- [`SetFeaturesV5()`](https://mianaz.github.io/scConvert/reference/SetFeaturesV5.md)
  : Set V5 feature index in HDF5 group

- [`SeuratLayerToAnnData()`](https://mianaz.github.io/scConvert/reference/SeuratLayerToAnnData.md)
  : Map Seurat layer name to AnnData layer name

- [`IndexToPointer()`](https://mianaz.github.io/scConvert/reference/SparsePointers.md)
  [`PointerToIndex()`](https://mianaz.github.io/scConvert/reference/SparsePointers.md)
  : Convert sparse matrix pointers to indices and vice versa

- [`SparseWrite()`](https://mianaz.github.io/scConvert/reference/SparseWrite.md)
  : Write a sparse matrix to an HDF5 dataset

- [`StringType()`](https://mianaz.github.io/scConvert/reference/StringType.md)
  : Generate an HDF5 string dtype

- [`IsDataFrame(`*`<H5D>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  [`IsDataFrame(`*`<H5Group>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  [`IsFactor(`*`<H5D>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  [`IsFactor(`*`<H5Group>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  [`IsList(`*`<H5Group>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  [`IsLogical(`*`<H5D>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  [`IsMatrix(`*`<H5D>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  [`IsMatrix(`*`<H5Group>`*`)`](https://mianaz.github.io/scConvert/reference/TestH5.md)
  : Test HDF5 datasets and groups to see what kind of data they are

- [`IsDataFrame()`](https://mianaz.github.io/scConvert/reference/TestObject.md)
  [`IsFactor()`](https://mianaz.github.io/scConvert/reference/TestObject.md)
  [`IsList()`](https://mianaz.github.io/scConvert/reference/TestObject.md)
  [`IsLogical()`](https://mianaz.github.io/scConvert/reference/TestObject.md)
  [`IsMatrix()`](https://mianaz.github.io/scConvert/reference/TestObject.md)
  : Test an object's class

- [`FormatTime()`](https://mianaz.github.io/scConvert/reference/Timestamp.md)
  [`Timestamp()`](https://mianaz.github.io/scConvert/reference/Timestamp.md)
  [`TSFormats()`](https://mianaz.github.io/scConvert/reference/Timestamp.md)
  : Create and work with timestamps

- [`TransferMetadataV5()`](https://mianaz.github.io/scConvert/reference/TransferMetadataV5.md)
  : Transfer metadata from V5 h5Seurat to h5ad

- [`UpdateKey()`](https://mianaz.github.io/scConvert/reference/UpdateKey.md)
  : Update a Seurat key

- [`UpdateSlots()`](https://mianaz.github.io/scConvert/reference/UpdateSlots.md)
  : Update slots in an object

- [`UtilsAssayCompat`](https://mianaz.github.io/scConvert/reference/UtilsAssayCompat.md)
  : Assay Compatibility Utility Functions

- [`UtilsH5Access`](https://mianaz.github.io/scConvert/reference/UtilsH5Access.md)
  : HDF5 Access Utility Functions

- [`V5LayerSupport`](https://mianaz.github.io/scConvert/reference/V5LayerSupport.md)
  : V5 Layer Support Functions

- [`LoomValidate0.1()`](https://mianaz.github.io/scConvert/reference/ValidateLoom.md)
  [`LoomValidate3.0.0()`](https://mianaz.github.io/scConvert/reference/ValidateLoom.md)
  : Validate Loom Files

- [`ValidateMultimodalObject()`](https://mianaz.github.io/scConvert/reference/ValidateMultimodalObject.md)
  : Validate Seurat object for multimodal export

- [`WriteH5SpatialV5()`](https://mianaz.github.io/scConvert/reference/WriteH5SpatialV5.md)
  : Handle spatial data for V5 objects

- [`WriteMode()`](https://mianaz.github.io/scConvert/reference/WriteMode.md)
  : Get the proper HDF5 connection mode for writing depending on
  overwrite status

- [`WriteV5Layer()`](https://mianaz.github.io/scConvert/reference/WriteV5Layer.md)
  : Write V5 Layer Data

- [`Writeable()`](https://mianaz.github.io/scConvert/reference/Writeable.md)
  : Is an HDF5 file or group writeable

- [`as.factor.H5D()`](https://mianaz.github.io/scConvert/reference/as.factor.H5D.md)
  : Convert an HDF5 dataset to a factor

- [`.add_spatialdata_coords()`](https://mianaz.github.io/scConvert/reference/dot-add_spatialdata_coords.md)
  : Add spatial coordinates from SpatialData to Seurat object

- [`.assay_to_modality()`](https://mianaz.github.io/scConvert/reference/dot-assay_to_modality.md)
  : Map Seurat assay name to h5mu modality name

- [`.attach_spatialdata_image()`](https://mianaz.github.io/scConvert/reference/dot-attach_spatialdata_image.md)
  : Attach an image array to a Seurat object

- [`.decode_vlen_utf8()`](https://mianaz.github.io/scConvert/reference/dot-decode_vlen_utf8.md)
  : Decode vlen-utf8 encoded raw data to character vector

- [`.encode_vlen_utf8()`](https://mianaz.github.io/scConvert/reference/dot-encode_vlen_utf8.md)
  : Encode character vector to vlen-utf8 raw data

- [`.fast_create_seurat()`](https://mianaz.github.io/scConvert/reference/dot-fast_create_seurat.md)
  : Fast Seurat object constructor (bypass CreateSeuratObject
  validation)

- [`.h5ad_write_dict_group()`](https://mianaz.github.io/scConvert/reference/dot-h5ad_write_dict_group.md)
  : Write an AnnData dict-style empty group

- [`.h5seurat_write_scaffold()`](https://mianaz.github.io/scConvert/reference/dot-h5seurat_write_scaffold.md)
  : Write required h5seurat root attributes and empty groups

- [`.modality_to_assay()`](https://mianaz.github.io/scConvert/reference/dot-modality_to_assay.md)
  : Map h5mu modality name to Seurat assay name

- [`.ome_ngff_to_hwc()`](https://mianaz.github.io/scConvert/reference/dot-ome_ngff_to_hwc.md)
  : Convert raw OME-NGFF array to height x width x channels format

- [`.raw_to_r_type()`](https://mianaz.github.io/scConvert/reference/dot-raw_to_r_type.md)
  : Convert raw bytes to R type

- [`.read_ome_ngff_image()`](https://mianaz.github.io/scConvert/reference/dot-read_ome_ngff_image.md)
  : Read an OME-NGFF image from a zarr directory

- [`.read_spatialdata_points()`](https://mianaz.github.io/scConvert/reference/dot-read_spatialdata_points.md)
  : Read points from a SpatialData points directory

- [`.read_spatialdata_shapes()`](https://mianaz.github.io/scConvert/reference/dot-read_spatialdata_shapes.md)
  : Read shapes from a SpatialData shapes directory

- [`.sd_read_column()`](https://mianaz.github.io/scConvert/reference/dot-sd_read_column.md)
  : Read a single column from a zarr DataFrame

- [`.soma_has_write_soma()`](https://mianaz.github.io/scConvert/reference/dot-soma_has_write_soma.md)
  : Check whether tiledbsoma::write_soma is available

- [`.soma_read_extra_measurements()`](https://mianaz.github.io/scConvert/reference/dot-soma_read_extra_measurements.md)
  : Read additional measurements as extra assays

- [`.soma_read_obs()`](https://mianaz.github.io/scConvert/reference/dot-soma_read_obs.md)
  : Read obs from a SOMA experiment

- [`.soma_read_obsm()`](https://mianaz.github.io/scConvert/reference/dot-soma_read_obsm.md)
  : Read obsm embeddings from a SOMA measurement

- [`.soma_read_obsp()`](https://mianaz.github.io/scConvert/reference/dot-soma_read_obsp.md)
  : Read obsp graphs from a SOMA measurement

- [`.soma_read_var()`](https://mianaz.github.io/scConvert/reference/dot-soma_read_var.md)
  : Read var from a SOMA measurement

- [`.soma_read_x()`](https://mianaz.github.io/scConvert/reference/dot-soma_read_x.md)
  : Read X matrix from a SOMA measurement

- [`.stream_composite_via_h5seurat()`](https://mianaz.github.io/scConvert/reference/dot-stream_composite_via_h5seurat.md)
  : Run two-step streaming conversion via a temp h5seurat file

- [`.stream_dense_matrix()`](https://mianaz.github.io/scConvert/reference/dot-stream_dense_matrix.md)
  : Stream a dense matrix from one HDF5 file to another

- [`.stream_h5ad_df_to_h5seurat_md()`](https://mianaz.github.io/scConvert/reference/dot-stream_h5ad_df_to_h5seurat_md.md)
  : Stream an h5ad-style DataFrame group to h5seurat meta.data

- [`.stream_h5mu_to_h5seurat()`](https://mianaz.github.io/scConvert/reference/dot-stream_h5mu_to_h5seurat.md)
  : Stream h5mu to h5seurat (direct HDF5-to-HDF5)

- [`.stream_h5seurat_md_to_h5ad_df()`](https://mianaz.github.io/scConvert/reference/dot-stream_h5seurat_md_to_h5ad_df.md)
  : Stream h5seurat meta.data to an h5ad-style DataFrame group

- [`.stream_h5seurat_to_h5mu()`](https://mianaz.github.io/scConvert/reference/dot-stream_h5seurat_to_h5mu.md)
  : Stream h5seurat to h5mu (direct HDF5-to-HDF5)

- [`.stream_h5seurat_to_loom()`](https://mianaz.github.io/scConvert/reference/dot-stream_h5seurat_to_loom.md)
  : Stream h5seurat to loom (direct HDF5-to-HDF5)

- [`.stream_loom_to_h5seurat()`](https://mianaz.github.io/scConvert/reference/dot-stream_loom_to_h5seurat.md)
  : Stream loom to h5seurat (direct HDF5-to-HDF5)

- [`.stream_sparse_group()`](https://mianaz.github.io/scConvert/reference/dot-stream_sparse_group.md)
  : Stream a sparse matrix group from one HDF5 file to another

- [`.writeH5AD_c()`](https://mianaz.github.io/scConvert/reference/dot-writeH5AD_c.md)
  : Fast C-based h5ad writer

- [`.writeH5Seurat_c()`](https://mianaz.github.io/scConvert/reference/dot-writeH5Seurat_c.md)
  : Fast C-based h5Seurat writer

- [`.write_ome_ngff_image()`](https://mianaz.github.io/scConvert/reference/dot-write_ome_ngff_image.md)
  : Write an image as an OME-NGFF zarr store

- [`.write_spatialdata_points()`](https://mianaz.github.io/scConvert/reference/dot-write_spatialdata_points.md)
  : Write point coordinates as a SpatialData points directory

- [`.write_spatialdata_shapes()`](https://mianaz.github.io/scConvert/reference/dot-write_spatialdata_shapes.md)
  : Write spatial coordinates as a SpatialData shapes directory

- [`.zarr_chunk_path_v2()`](https://mianaz.github.io/scConvert/reference/dot-zarr_chunk_path_v2.md)
  : Build chunk file path for zarr v2

- [`.zarr_compress()`](https://mianaz.github.io/scConvert/reference/dot-zarr_compress.md)
  : Compress chunk data

- [`.zarr_create_group()`](https://mianaz.github.io/scConvert/reference/dot-zarr_create_group.md)
  : Create a zarr v2 group directory with metadata

- [`.zarr_decompress()`](https://mianaz.github.io/scConvert/reference/dot-zarr_decompress.md)
  : Decompress chunk data

- [`.zarr_get_n_obs()`](https://mianaz.github.io/scConvert/reference/dot-zarr_get_n_obs.md)
  : Get number of observations from zarr store

- [`.zarr_get_n_vars()`](https://mianaz.github.io/scConvert/reference/dot-zarr_get_n_vars.md)
  : Get number of variables from zarr store

- [`.zarr_list_children()`](https://mianaz.github.io/scConvert/reference/dot-zarr_list_children.md)
  : List children of a zarr group

- [`.zarr_node_type()`](https://mianaz.github.io/scConvert/reference/dot-zarr_node_type.md)
  : Determine zarr node type

- [`.zarr_parallel_write()`](https://mianaz.github.io/scConvert/reference/dot-zarr_parallel_write.md)
  : Write multiple zarr arrays in parallel

- [`.zarr_parse_dtype()`](https://mianaz.github.io/scConvert/reference/dot-zarr_parse_dtype.md)
  : Parse zarr v2 dtype string to R type info

- [`.zarr_read_anndata_column()`](https://mianaz.github.io/scConvert/reference/dot-zarr_read_anndata_column.md)
  : Read an AnnData DataFrame column from zarr store

- [`.zarr_read_anndata_index()`](https://mianaz.github.io/scConvert/reference/dot-zarr_read_anndata_index.md)
  : Read AnnData index from obs or var group

- [`.zarr_read_anndata_index_from_attrs()`](https://mianaz.github.io/scConvert/reference/dot-zarr_read_anndata_index_from_attrs.md)
  : Read an anndata index given pre-read attrs

- [`.zarr_read_anndata_matrix()`](https://mianaz.github.io/scConvert/reference/dot-zarr_read_anndata_matrix.md)
  : Read an AnnData matrix (dense or sparse) from zarr store

- [`.zarr_read_attrs()`](https://mianaz.github.io/scConvert/reference/dot-zarr_read_attrs.md)
  : Read zarr node attributes

- [`.zarr_read_chunked()`](https://mianaz.github.io/scConvert/reference/dot-zarr_read_chunked.md)
  : Read multi-chunk zarr v2 numeric array

- [`.zarr_read_json()`](https://mianaz.github.io/scConvert/reference/dot-zarr_read_json.md)
  : Read and parse a JSON file

- [`.zarr_read_numeric()`](https://mianaz.github.io/scConvert/reference/dot-zarr_read_numeric.md)
  : Read a zarr v2 numeric array from disk

- [`.zarr_read_strings()`](https://mianaz.github.io/scConvert/reference/dot-zarr_read_strings.md)
  : Read a vlen-utf8 encoded string array from zarr v2 store

- [`.zarr_store_version()`](https://mianaz.github.io/scConvert/reference/dot-zarr_store_version.md)
  : Detect zarr store format version

- [`.zarr_write_anndata_dataframe()`](https://mianaz.github.io/scConvert/reference/dot-zarr_write_anndata_dataframe.md)
  : Write an AnnData DataFrame (obs or var) to zarr store

- [`.zarr_write_attrs()`](https://mianaz.github.io/scConvert/reference/dot-zarr_write_attrs.md)
  : Write zarr attributes to .zattrs file

- [`.zarr_write_numeric()`](https://mianaz.github.io/scConvert/reference/dot-zarr_write_numeric.md)
  : Write a zarr v2 numeric array

- [`.zarr_write_strings()`](https://mianaz.github.io/scConvert/reference/dot-zarr_write_strings.md)
  : Write a zarr v2 vlen-utf8 string array

- [`DefaultAssay(`*`<h5SI>`*`)`](https://mianaz.github.io/scConvert/reference/h5SI.md)
  [`print(`*`<h5SI>`*`)`](https://mianaz.github.io/scConvert/reference/h5SI.md)
  : Tools for handling h5Seurat indexes

- [`scAppendData()`](https://mianaz.github.io/scConvert/reference/scAppendData.md)
  :

  Append data from an h5Seurat file to a preexisting
  [`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html)
  object

- [`sc_find_cli()`](https://mianaz.github.io/scConvert/reference/sc_find_cli.md)
  : Find the scconvert CLI binary
