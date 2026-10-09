use core::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BadConversion;

pub fn convert_substring<T>(value: &str, i: usize, j: usize) -> Result<T, BadConversion>
where
    T: FromStr,
{
    let slice = substring_1_based(value, i, j)?;
    let trimmed = slice.trim_start();
    trimmed.parse::<T>().map_err(|_| BadConversion)
}

pub fn convert_substring_unsigned(value: &str, i: usize, j: usize) -> Result<u32, BadConversion> {
    let tmp = convert_substring::<i32>(value, i, j)?;
    if tmp < 0 {
        return Err(BadConversion);
    }
    Ok(tmp as u32)
}

pub fn substring_is_blank(value: &str, i: usize, j: usize) -> Result<bool, BadConversion> {
    let slice = substring_1_based(value, i, j)?;
    Ok(slice.chars().all(char::is_whitespace))
}

fn substring_1_based(value: &str, i: usize, j: usize) -> Result<&str, BadConversion> {
    if i < 1
        || i > j + 1
        || j > value.len()
        || !value.is_char_boundary(i - 1)
        || !value.is_char_boundary(j)
    {
        return Err(BadConversion);
    }
    Ok(&value[(i - 1)..j])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_one_based_fixed_columns() {
        let value = "ATOM    42  ";
        assert_eq!(convert_substring::<i32>(value, 7, 10).unwrap(), 42);
        assert!(substring_is_blank(value, 11, 12).unwrap());
    }

    #[test]
    fn rejects_negative_unsigned_like_reference_specialization() {
        assert!(convert_substring_unsigned(" -3", 1, 3).is_err());
        assert_eq!(convert_substring_unsigned("  3", 1, 3).unwrap(), 3);
    }
}
