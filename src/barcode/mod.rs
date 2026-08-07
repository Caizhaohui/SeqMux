pub mod config;
pub mod matcher;

pub use config::{
    load_barcodes_csv, parse_barcodes_csv, sample_key_for, BarcodeConfig, CompiledFivePrime,
    CompiledThreePrime, DemuxMode, SampleKey, TsoPattern,
};
pub use matcher::{
    extract_five_umi, extract_three_umi, hamming, match_five_prime, match_three_prime, MatchResult,
};
