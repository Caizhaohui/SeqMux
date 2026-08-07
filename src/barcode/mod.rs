pub mod config;
pub mod matcher;

pub use config::{
    load_barcodes_csv, parse_barcodes_csv, BarcodeConfig, CompiledBarcode, DemuxMode, SampleEntry,
    SampleKey, TsoPattern,
};
pub use matcher::{
    assign_paired, assign_single, barcode_distance_3p, barcode_distance_5p, extract_umi_3p,
    extract_umi_5p, hamming, sample_key, MatchResult,
};
