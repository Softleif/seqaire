//! Inspect SAM flags via the [`BamFlags`] newtype. Named constants (e.g. [`FLAG_UNMAPPED`],
//! [`FLAG_REVERSE`]) prevent accidental use of raw `u16` literals.

// r[impl io.named_constants]

// r[impl io.typed_flags]
/// Strongly-typed BAM flags that prevent accidentally mixing raw u16 values.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct BamFlags(u16);

impl BamFlags {
    pub const fn new(raw: u16) -> Self {
        Self(raw)
    }

    pub const fn raw(self) -> u16 {
        self.0
    }

    pub const fn is_set(self, flag: u16) -> bool {
        self.0 & flag != 0
    }

    pub const fn is_paired(self) -> bool {
        self.is_set(FLAG_PAIRED)
    }
    pub const fn is_proper_pair(self) -> bool {
        self.is_set(FLAG_PROPER_PAIR)
    }
    pub const fn is_unmapped(self) -> bool {
        self.is_set(FLAG_UNMAPPED)
    }
    pub const fn is_reverse(self) -> bool {
        self.is_set(FLAG_REVERSE)
    }
    pub const fn is_first_in_template(self) -> bool {
        self.is_set(FLAG_FIRST_IN_TEMPLATE)
    }
    pub const fn is_second_in_template(self) -> bool {
        self.is_set(FLAG_SECOND_IN_TEMPLATE)
    }
    pub const fn is_secondary(self) -> bool {
        self.is_set(FLAG_SECONDARY)
    }
    pub const fn is_duplicate(self) -> bool {
        self.is_set(FLAG_DUPLICATE)
    }
    pub const fn is_supplementary(self) -> bool {
        self.is_set(FLAG_SUPPLEMENTARY)
    }
}

/// Read is paired in sequencing.
pub const FLAG_PAIRED: u16 = 0x1;
/// Read is mapped in a proper pair.
pub const FLAG_PROPER_PAIR: u16 = 0x2;
/// Read is unmapped.
pub const FLAG_UNMAPPED: u16 = 0x4;
/// Mate is unmapped.
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
/// Read is mapped to the reverse strand.
pub const FLAG_REVERSE: u16 = 0x10;
/// Mate is mapped to the reverse strand.
pub const FLAG_MATE_REVERSE: u16 = 0x20;
/// Read is the first in a template (pair).
pub const FLAG_FIRST_IN_TEMPLATE: u16 = 0x40;
/// Read is the second (last) in a template (pair).
pub const FLAG_SECOND_IN_TEMPLATE: u16 = 0x80;
/// Alignment is secondary (not primary).
pub const FLAG_SECONDARY: u16 = 0x100;
/// Read fails platform/vendor quality checks.
pub const FLAG_FAILED_QC: u16 = 0x200;
/// Read is a PCR or optical duplicate.
pub const FLAG_DUPLICATE: u16 = 0x400;
/// Alignment is supplementary.
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;
