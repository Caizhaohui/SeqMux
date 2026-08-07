/// Reverse-complement a DNA sequence. Non-ACGT bases become `N`.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        out.push(match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        });
    }
    out
}

/// Normalize a FASTQ read ID for pair matching: strip `/1` `/2` and
/// Illumina space-separated second field (` read:N:N:...`).
pub fn normalize_read_id(name: &[u8]) -> &[u8] {
    let mut end = name.len();
    // Strip after first space or tab
    if let Some(pos) = name.iter().position(|&b| b == b' ' || b == b'\t') {
        end = pos;
    }
    let core = &name[..end];
    if core.len() >= 2 {
        let last2 = &core[core.len() - 2..];
        if last2 == b"/1" || last2 == b"/2" {
            return &core[..core.len() - 2];
        }
    }
    core
}

/// Sanitize a string for use as a filename component.
pub fn sanitize_filename_component(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        if c.is_ascii_alphanumeric() || c == '-' || c == '_' || c == '.' {
            out.push(c);
        } else {
            out.push('_');
        }
    }
    if out.is_empty() {
        out.push_str("unnamed");
    }
    // Windows reserved names
    let upper = out.to_ascii_uppercase();
    match upper.as_str() {
        "CON" | "PRN" | "AUX" | "NUL" | "COM1" | "COM2" | "COM3" | "COM4" | "COM5" | "COM6"
        | "COM7" | "COM8" | "COM9" | "LPT1" | "LPT2" | "LPT3" | "LPT4" | "LPT5" | "LPT6"
        | "LPT7" | "LPT8" | "LPT9" => {
            format!("_{out}")
        }
        _ => out,
    }
}

/// Format a large integer with thousands separators.
pub fn format_count(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out.chars().rev().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_revcomp() {
        assert_eq!(revcomp(b"ACGTN"), b"NACGT");
        assert_eq!(revcomp(b"ACGT"), b"ACGT");
        assert_eq!(revcomp(b"AAAA"), b"TTTT");
        assert_eq!(revcomp(b"acgt"), b"ACGT");
    }

    #[test]
    fn test_revcomp_idempotent() {
        let seq = b"ATCGGATCNNN";
        assert_eq!(revcomp(&revcomp(seq)), seq.to_vec());
    }

    #[test]
    fn test_normalize_read_id() {
        assert_eq!(normalize_read_id(b"read1/1"), b"read1");
        assert_eq!(normalize_read_id(b"read1/2"), b"read1");
        assert_eq!(
            normalize_read_id(b"A01234:87:HXXX:1:1101:1000:1000 1:N:0:ATC"),
            b"A01234:87:HXXX:1:1101:1000:1000"
        );
    }

    #[test]
    fn test_sanitize_filename() {
        assert_eq!(sanitize_filename_component("sample-1"), "sample-1");
        assert_eq!(sanitize_filename_component("a/b:c"), "a_b_c");
        assert_eq!(sanitize_filename_component("CON"), "_CON");
    }
}
