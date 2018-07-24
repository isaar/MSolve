namespace ISAAR.MSolve.LinearAlgebra.Output.Formatting
{
    /// <summary>
    /// For each number, a fixed number of decimal points are displayed, e.g -1234.56. Also controls alignment by reserving 
    /// enough space for the max number of integer digits expected.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class FixedPointFormat: INumericFormat
    {
        /// <summary>
        /// Decimal digits to display, not counting the decimal separator (. or ,).
        /// </summary>
        public int NumDecimalDigits { get; set; } = 6;

        /// <summary>
        /// Not counting the sign (+ or -) or the decimal separator (. or ,). This property is used to reserve enough space, so 
        /// that matrix/vector entries in the same column can be aligned. If the <see cref="MaxIntegerDigits"/> chosen is not 
        /// enough, then the alignment will be destroyed, but the number will be displayed in full. If this is unwanted and the 
        /// maximum number of integer digits of the entries of the matrix/vector cannot be guessed, consider using 
        /// <see cref="ExponentialFormat"/> instead.
        /// </summary>
        public int MaxIntegerDigits { get; set; } = 4;

        /// <summary>
        /// Returns a tring that will be used as a template for formatting numbers during output.
        /// </summary>
        public string GetRealNumberFormat()
        {
            int totalWidth = NumDecimalDigits + MaxIntegerDigits + 2; // +1 for sign, + 1 for comma/period
            return "{0," + $"{totalWidth}:F{NumDecimalDigits}" + "}";
        }
    }
}
