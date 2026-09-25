pub mod adapter;
pub mod quality;

pub use adapter::{
    find_3p_adapter, find_3p_adapter_fast, find_3p_adapter_generic, trim_3p_adapter, AdapterMatch,
    CompiledAdapter,
};
pub use quality::{nextseq_trim_index, quality_trim_bounds};
