use color_eyre::eyre::{Context, Report};
use std::{ops::Deref, str::FromStr};

/// A probability value guaranteed to be in the range [0, 1]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[must_use]
pub struct Probability(f64);

#[derive(Debug, thiserror::Error)]
#[error("Probability value `{0}` is outside the valid range [0, 1]")]
pub struct ProbabilityError(f64);

impl Probability {
    /// Creates a new Probability, returning an error if the value is outside [0, 1]
    pub const fn new(value: f64) -> Result<Self, ProbabilityError> {
        if value < 0.0 || value > 1.0 || !value.is_finite() {
            Err(ProbabilityError(value))
        } else {
            Ok(Self(value))
        }
    }

    /// Creates a new Probability, panicking if the value is outside [0, 1]
    ///
    /// Can be used to build Probability in const contexts.
    #[expect(clippy::panic, reason = "panicking is intentional here")]
    pub const fn new_panicky(value: f64) -> Self {
        match Self::new(value) {
            Ok(p) => p,
            Err(_e) => panic!("invalid probability"),
        }
    }

    /// Returns the inverted probability (1 - p)
    pub const fn inverted(self) -> Self {
        Self(1.0 - self.0)
    }

    /// The probability value 0.0
    pub const ZERO: Self = Self(0.0);
    /// The probability value 1.0
    pub const ONE: Self = Self(1.0);
}

impl std::fmt::Display for Probability {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(precision) = f.precision() {
            // If we received a precision, we use it.
            write!(f, "{1:.*}", precision, self.0)
        } else {
            // Otherwise we default to 2.
            write!(f, "{:.2}", self.0)
        }
    }
}

impl FromStr for Probability {
    type Err = Report;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.parse::<f64>()
            .wrap_err("Failed to parse probability value")
            .and_then(|v| Self::new(v).wrap_err("Probability value must be in the range [0, 1]"))
            .wrap_err_with(|| format!("Invalid probability: {s}"))
    }
}

impl Deref for Probability {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl TryFrom<f64> for Probability {
    type Error = ProbabilityError;

    fn try_from(value: f64) -> Result<Self, Self::Error> {
        Self::new(value)
    }
}

impl Default for Probability {
    fn default() -> Self {
        Self(0.0)
    }
}

impl serde::Serialize for Probability {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.0.serialize(serializer)
    }
}

impl<'de> serde::Deserialize<'de> for Probability {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let value = f64::deserialize(deserializer)?;
        Self::new(value).map_err(serde::de::Error::custom)
    }
}

#[test]
fn fails_outside_range() {
    assert!(Probability::new(-0.1).is_err());
    assert!(Probability::new(1.1).is_err());
    assert!(Probability::new(f64::NAN).is_err());
    assert!(Probability::new(f64::INFINITY).is_err());
    assert!(Probability::new(f64::NEG_INFINITY).is_err());
}
