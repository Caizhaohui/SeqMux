pub mod config;
pub mod matcher;

pub use config::{
    load_barcodes_csv, pack_8bp_2bit, parse_barcodes_csv, BarcodeConfig, CompiledBarcode,
    DemuxMode, ExactDual8Matcher, OrientationMode, SampleEntry, SampleKey, TsoPattern,
};
pub use matcher::{
    assign_paired, assign_single, barcode_distance_3p, barcode_distance_5p, extract_umi_3p,
    extract_umi_5p, hamming, sample_key, MatchResult,
};
