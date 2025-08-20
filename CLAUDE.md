# Ptera Software Development Guidelines for Claude

## Writing Style Guidelines

### Terminology
- **"Ptera Software"**: When writing as text, always write as two words without a 
hyphen, each being capitalized (never "ptera", "ptera software", or "PteraSoftware"). 
When writing as a package to be imported, use ```import pterasoftware as ps```  
- **"cross section"**: Always write as two words, never hyphenated 
(not "cross-section")  
- **Object references**: When referring to code objects, use proper class naming 
convention. This will illustrate that we are talking about a code object, not an 
abstraction. Therefore, don't include the "object" suffix (e.g. "this WingCrossSection 
object") unless it is needed for clarity. In summary, when talking about code objects:
  - ✅ "the previous WingCrossSection"
  - ❌ "the previous cross section"
  - ✅ "this Wing"
  - ❌ "this wing"
- **Abstract references**: When referring to abstractions, use lowercase and separate 
individual words with a space (e.g. "an airplane's wings are used to generate lift" and 
"the cross section of a wing typically has a streamlined shape known as an airfoil"). 
This is to distinguish them from code objects.

### Docstring Style
- Follow existing PteraSoftware docstring conventions  
- Use reStructuredText (rST) formatting guidelines  
- Maintain consistent parameter descriptions  
- Preserve existing documentation structure and completeness unless we are explicitly 
updating or improving it  
- Always include units in parameter descriptions where applicable  
- Keep parameter descriptions complete.  
- If a parameter is a numpy array, specify the expected shape and data type. For  
example, say "(3,) ndarray of floats" for a 1D array of 3 floats.

### Code Formatting
- Follow existing code style (black) and conventions  
- Maintain consistent indentation and spacing  
- Preserve existing comment structure and detail level  

### Variable Naming
- Use descriptive variables names that clearly indicate their purpose  
- Use lowercase with underscores for variable names
- Variables that are coordinates should be named 1D ndarrays, and have a suffix 
abbreviation indicating their reference frame. The reference frames used are:  
  - The wind frame (wind_frame)
  - The geometry frame (geometry_frame)
  - The wing frame (wing_frame)
  - The wing cross section frame (wing_cross_section_frame)
  - The airfoil frame (airfoil_frame)
  - The for a 1D ndarray of x-coordinates in the body frame, or "x_wing" for a 1D ndarray of x-coordinates in the wing frame.

### Miscellaneous Guidelines
- Use clear, descriptive variable names  
- Avoid abbreviations unless they are well-known in the context  
- In docstrings and comments, never use em-dashes (—) or en-dashes (–); always use 
hyphens (-) for clarity  
- Never use emojis in code, comments, or docstrings  