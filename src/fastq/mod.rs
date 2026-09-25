pub mod reader;
pub mod record;
pub mod writer;

pub use reader::{
    validate_fastq, ConcurrentPairedReader, FastqReader, FastqSource, InputReader,
    PairedFastqReader, ValidateStats,
};
pub use record::{OwnedFastqRecord, ReadOrPair, ReadPair};
pub use writer::{check_output_path, FastqWriter};
