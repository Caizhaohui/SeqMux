pub mod adapter;
pub mod quality;

pub use adapter::{find_3p_adapter, trim_3p_adapter, AdapterMatch};
pub use quality::{nextseq_trim_index, quality_trim_bounds};
