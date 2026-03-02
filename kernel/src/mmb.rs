// mmb.rs - Complete MMB file reader, parser, and inspector
// Use as:
//   - Library: mod mmb; use mmb::MmbFile;
//   - Binary: cargo run --bin mmb -- [--full] <file.mmb>

use std::env;
use std::fs;
use std::io;

// ============================================================================
// LIBRARY PART - Use in your code
// ============================================================================

/// MMB file header structure (40 bytes)
#[derive(Debug, Clone)]
pub struct MmbHeader {
    pub magic: [u8; 4],
    pub version: u8,
    pub num_sorts: u8,
    pub reserved: u16,
    pub num_terms: u32,
    pub num_thms: u32,
    pub p_terms: u32,
    pub p_thms: u32,
    pub p_proof: u32,
    pub reserved2: u32,
    pub p_index: u64,
}

impl MmbHeader {
    pub fn from_bytes(buf: &[u8]) -> Result<Self, &'static str> {
        if buf.len() < 40 {
            return Err("Buffer too small for header");
        }

        Ok(MmbHeader {
            magic: [buf[0], buf[1], buf[2], buf[3]],
            version: buf[4],
            num_sorts: buf[5],
            reserved: u16::from_le_bytes([buf[6], buf[7]]),
            num_terms: u32::from_le_bytes([buf[8], buf[9], buf[10], buf[11]]),
            num_thms: u32::from_le_bytes([buf[12], buf[13], buf[14], buf[15]]),
            p_terms: u32::from_le_bytes([buf[16], buf[17], buf[18], buf[19]]),
            p_thms: u32::from_le_bytes([buf[20], buf[21], buf[22], buf[23]]),
            p_proof: u32::from_le_bytes([buf[24], buf[25], buf[26], buf[27]]),
            reserved2: u32::from_le_bytes([buf[28], buf[29], buf[30], buf[31]]),
            p_index: u64::from_le_bytes([
                buf[32], buf[33], buf[34], buf[35], buf[36], buf[37], buf[38], buf[39],
            ]),
        })
    }

    pub fn validate(&self) -> Result<(), String> {
        if &self.magic != b"MM0B" {
            return Err(format!(
                "Invalid magic: expected 'MM0B', got {:?}",
                self.magic
            ));
        }
        if self.version != 1 {
            return Err(format!("Unsupported version: {}", self.version));
        }
        if self.p_terms >= self.p_thms || self.p_thms >= self.p_proof {
            return Err("Invalid pointer ordering".to_string());
        }
        Ok(())
    }
}

/// Sort data with modifiers
#[derive(Debug, Clone)]
pub struct SortData {
    pub raw: u8,
    pub pure: bool,
    pub strict: bool,
    pub provable: bool,
    pub free: bool,
}

impl SortData {
    pub fn from_byte(byte: u8) -> Self {
        SortData {
            raw: byte,
            pure: (byte & 0x01) != 0,
            strict: (byte & 0x02) != 0,
            provable: (byte & 0x04) != 0,
            free: (byte & 0x08) != 0,
        }
    }
}

/// Main MMB file structure
pub struct MmbFile {
    pub header: MmbHeader,
    pub data: Vec<u8>,
    pub sorts: Vec<SortData>,
}

impl MmbFile {
    pub fn load(filename: &str) -> io::Result<Self> {
        let data = fs::read(filename)?;
        Self::from_bytes(data)
    }

    pub fn from_bytes(data: Vec<u8>) -> io::Result<Self> {
        let header = MmbHeader::from_bytes(&data)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        header
            .validate()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let mut sorts = Vec::new();
        let sorts_start = 40;
        for i in 0..header.num_sorts as usize {
            if sorts_start + i < data.len() {
                sorts.push(SortData::from_byte(data[sorts_start + i]));
            }
        }

        Ok(MmbFile {
            header,
            data,
            sorts,
        })
    }

    pub fn proof_stream(&self) -> &[u8] {
        let start = self.header.p_proof as usize;
        let end = if self.header.p_index > 0 {
            self.header.p_index as usize
        } else {
            self.data.len()
        };
        &self.data[start..end]
    }

    // Quick summary
    pub fn print_summary(&self) {
        println!("=== MMB File Summary ===");
        println!(
            "Magic:    {:?}",
            std::str::from_utf8(&self.header.magic).unwrap_or("?")
        );
        println!("Version:  {}", self.header.version);
        println!("Sorts:    {}", self.header.num_sorts);
        println!("Terms:    {}", self.header.num_terms);
        println!("Theorems: {}", self.header.num_thms);
        println!("Size:     {} bytes", self.data.len());
    }

    // Detailed inspection
    pub fn print_detailed(&self) {
        println!("=== COMPLETE MMB FILE DUMP ===");
        println!("Size: {} bytes (0x{:x})", self.data.len(), self.data.len());
        println!();

        // Header
        println!("=== HEADER (0x00-0x27) ===");
        println!(
            "0x00-03 Magic:     {:?} \"{}\"",
            self.header.magic,
            std::str::from_utf8(&self.header.magic).unwrap_or("?")
        );
        println!("0x04    Version:   {}", self.header.version);
        println!("0x05    Sorts:     {}", self.header.num_sorts);
        println!("0x06-07 Reserved:  {:#06x}", self.header.reserved);
        println!("0x08-0b #Terms:    {}", self.header.num_terms);
        println!("0x0c-0f #Thms:     {}", self.header.num_thms);
        println!(
            "0x10-13 p_terms:   {:#010x} (byte {})",
            self.header.p_terms, self.header.p_terms
        );
        println!(
            "0x14-17 p_thms:    {:#010x} (byte {})",
            self.header.p_thms, self.header.p_thms
        );
        println!(
            "0x18-1b p_proof:   {:#010x} (byte {})",
            self.header.p_proof, self.header.p_proof
        );
        println!("0x1c-1f Reserved2: {:#010x}", self.header.reserved2);
        println!(
            "0x20-27 p_index:   {:#018x} (byte {})",
            self.header.p_index, self.header.p_index
        );

        // Validation
        println!("\n=== VALIDATION ===");
        if &self.header.magic == b"MM0B" {
            println!("✓ Valid magic number");
        } else {
            println!("❌ Invalid magic");
        }

        if self.header.version == 1 {
            println!("✓ Valid version");
        } else {
            println!("⚠ Unexpected version");
        }

        // Check alignment
        let align_to_8 = |n: usize| (n + 7) & !7;
        let expected_p_terms = align_to_8(40 + self.header.num_sorts as usize);
        if self.header.p_terms as usize == expected_p_terms {
            println!("✓ p_terms correctly aligned to 8-byte boundary");
        } else {
            println!(
                "⚠ p_terms = {}, expected {} (alignment)",
                self.header.p_terms, expected_p_terms
            );
        }

        if self.header.p_terms < self.header.p_thms && self.header.p_thms < self.header.p_proof {
            println!("✓ Pointers in correct order");
        }

        // Sorts
        println!(
            "\n=== SORTS ({} bytes at offset 40) ===",
            self.header.num_sorts
        );
        for (i, sort) in self.sorts.iter().enumerate() {
            print!("  Sort {:2}: {:#04x} ", i, sort.raw);
            let mut flags = Vec::new();
            if sort.pure {
                flags.push("pure");
            }
            if sort.strict {
                flags.push("strict");
            }
            if sort.provable {
                flags.push("provable");
            }
            if sort.free {
                flags.push("free");
            }
            println!(
                "{}",
                if flags.is_empty() {
                    "(none)".to_string()
                } else {
                    flags.join(" ")
                }
            );
        }

        // Padding
        let sorts_end = 40 + self.header.num_sorts as usize;
        if sorts_end < self.header.p_terms as usize {
            let padding_size = self.header.p_terms as usize - sorts_end;
            println!(
                "\n=== PADDING (0x{:x}-0x{:x}) ===",
                sorts_end,
                self.header.p_terms - 1
            );
            println!(
                "Alignment padding: {} bytes (for 8-byte alignment)",
                padding_size
            );
        }

        // Terms
        println!(
            "\n=== TERMS TABLE (0x{:x}-0x{:x}) ===",
            self.header.p_terms,
            self.header.p_thms - 1
        );
        println!("Number of entries: {}", self.header.num_terms);

        let max_show = self.header.num_terms.min(10);
        if max_show > 0 {
            println!("\nFirst {} term entries:", max_show);
        }
        for i in 0..max_show {
            let offset = self.header.p_terms as usize + (i as usize * 8);
            if offset + 8 <= self.data.len() {
                let entry = &self.data[offset..offset + 8];
                let num_args = entry[0];
                let sort = entry[1];
                let p_args = u32::from_le_bytes([entry[4], entry[5], entry[6], entry[7]]);
                println!(
                    "  Term {}: {:02x?}  (sort={}, args={}, p_args=0x{:x})",
                    i, entry, sort, num_args, p_args
                );
            }
        }
        if self.header.num_terms > max_show {
            println!("  ... {} more terms", self.header.num_terms - max_show);
        }

        // Theorems
        println!(
            "\n=== THEOREMS TABLE (0x{:x}-0x{:x}) ===",
            self.header.p_thms,
            self.header.p_proof - 1
        );
        println!("Number of entries: {}", self.header.num_thms);

        let max_show = self.header.num_thms.min(10);
        if max_show > 0 {
            println!("\nFirst {} theorem entries:", max_show);
        }
        for i in 0..max_show {
            let offset = self.header.p_thms as usize + (i as usize * 8);
            if offset + 8 <= self.data.len() {
                let entry = &self.data[offset..offset + 8];
                let num_args = entry[0];
                let p_args = u32::from_le_bytes([entry[4], entry[5], entry[6], entry[7]]);
                println!(
                    "  Thm {}: {:02x?}  (args={}, p_args=0x{:x})",
                    i, entry, num_args, p_args
                );
            }
        }
        if self.header.num_thms > max_show {
            println!("  ... {} more theorems", self.header.num_thms - max_show);
        }

        // Proof stream
        let proof_end = if self.header.p_index > 0 {
            self.header.p_index as usize
        } else {
            self.data.len()
        };
        println!(
            "\n=== PROOF STREAM (0x{:x}-0x{:x}) ===",
            self.header.p_proof,
            proof_end - 1
        );
        println!("Size: {} bytes", proof_end - self.header.p_proof as usize);

        println!("\nFirst 64 bytes of proof stream:");
        self.hex_dump(self.header.p_proof as usize, 64);

        // Index
        if self.header.p_index > 0 {
            println!(
                "\n=== INDEX (0x{:x}-0x{:x}) ===",
                self.header.p_index,
                self.data.len() - 1
            );
            println!(
                "Size: {} bytes",
                self.data.len() - self.header.p_index as usize
            );
        }
    }

    pub fn hex_dump(&self, start: usize, len: usize) {
        let end = (start + len).min(self.data.len());

        for offset in (start..end).step_by(16) {
            let chunk_end = (offset + 16).min(end);
            let chunk = &self.data[offset..chunk_end];

            print!("{:08x}: ", offset);

            for (i, byte) in chunk.iter().enumerate() {
                print!("{:02x} ", byte);
                if i == 7 {
                    print!(" ");
                }
            }

            for i in chunk.len()..16 {
                print!("   ");
                if i == 7 {
                    print!(" ");
                }
            }

            print!(" |");
            for byte in chunk {
                let c = if *byte >= 0x20 && *byte <= 0x7e {
                    *byte as char
                } else {
                    '.'
                };
                print!("{}", c);
            }
            println!("|");
        }
    }
}

// ============================================================================
// BINARY PART - Command line tool
// ============================================================================

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} [--full] <file.mmb>", args[0]);
        eprintln!();
        eprintln!("Options:");
        eprintln!("  --full    Show complete detailed dump");
        eprintln!("  (default) Show summary only");
        std::process::exit(1);
    }

    let (full_mode, filename) = if args.len() == 3 && args[1] == "--full" {
        (true, &args[2])
    } else if args.len() == 2 {
        (false, &args[1])
    } else {
        eprintln!("Invalid arguments");
        std::process::exit(1);
    };

    let mmb = MmbFile::load(filename)?;

    if full_mode {
        mmb.print_detailed();
    } else {
        mmb.print_summary();
    }

    Ok(())
}
