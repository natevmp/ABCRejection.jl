# Gemini Instruction Profile: Array Naming Convention

## Core Naming Rules
1. **General Variables:** Use `camelCase`.
   - *Correct:* `refTime`, `totalValue`
   - *Incorrect:* `ref_time`, `Total_value`

2. **The Underscore Protocol (Arrays Only):**
   Underscores are strictly reserved for defining data-to-index relationships. This distinguishes flat arrays from nested arrays and specifies their layout.
   
   - **Indexing Array (Coordinate Vector):** Prefix the index name with a leading underscore.
     - *Example:* `_time`, `_t`, `_lat`
   - **Flat 1D Data Array:** `[element]_[index]`
     - *Example:* `pressure_time`, `p_t`
   - **Flat 2D Data Array:** `[element]_[index1]_[index2]` (all lowercase indices denote a flat contiguous array)
     - *Example:* `temp_x_y` (where `_x` and `_y` are the coordinate arrays)
   - **Flat ND Data Array:** `[element]_[index1]_[index2]_..._[indexN]`

   **Missing Index Arrays:** If a physical coordinate array does not exist, name the index `_id` or `_xid` (where `x` is the first letter of the entity, e.g., `_sid` for simulations). Casing rules are relaxed for index names.

   **Nested Arrays (Array of Arrays):**
   To distinguish nested structures (e.g., a vector of vectors) from flat multidimensional arrays:
   - A capital letter is introduced *only* when a **new level of encapsulation** is introduced.
   - The ordering of indices in the variable name goes from **innermost index to outermost index**.
   - *Example (1D nested in 1D):* A 1D array `element_bid` which is encapsulated inside an outer 1D array indexed by `_aid` is written as: `element_bid_Aid`.
   - *Example (1D nested in 2D flat, nested in 1D):* A 1D array `element_p` (indexed by `_p`) which is an element of a 2D flat array (indexed by `_a` and `_b`), which itself is an element of a 1D array (indexed by `_t`). Since the encapsulation boundaries happen at `a` and `t`, they are capitalized, resulting in: `element_p_A_b_T`.

## Operational Requirements
- **Strict Separation:** Never use underscores in standard variables. If a variable contains an underscore, it MUST represent an array following this protocol.
- **Do not rewrite existing code** unless explicitly asked to modify it.

## Validation Step
Before outputting code, verify:
1. Are flat 1D/ND arrays formatted with lowercase index names (e.g., `data_index`)?
2. Are encapsulation boundaries in nested arrays marked with a capital letter (e.g., `inner_Outer`)?
3. Are coordinate arrays prefixed with a leading underscore (e.g., `_t`)?

## Additional request
- Do not attempt to rewrite code that I have not specifically asked to be altered.