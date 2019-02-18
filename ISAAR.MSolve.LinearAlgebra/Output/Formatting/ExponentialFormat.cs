namespace ISAAR.MSolve.LinearAlgebra.Output.Formatting
{
    /// <summary>
    /// Numbers will be displayed in scientific (exponential) format: e.g. -1.05E+003. The caller can choose how many decimal 
    /// will be displayed and the resulting output will be aligned.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ExponentialFormat: INumericFormat
    {
        /// <summary>
        /// Decimal digits to display, not counting the decimal separator (. or ,).
        /// </summary>
        public int NumDecimalDigits { get; set; } = 6;

        /// <summary>
        /// Returns a tring that will be used as a template for formatting numbers during output.
        /// </summary>
        public string GetRealNumberFormat()
        {
            // -1.05E+003 => sign + integer + decimal separator + E + exponent sign + exponent digits + decimal digits.
            int totalWidth = NumDecimalDigits + 8;
            return "{0," + $"{totalWidth}:E{NumDecimalDigits}" + "}";
        }
    }
}
