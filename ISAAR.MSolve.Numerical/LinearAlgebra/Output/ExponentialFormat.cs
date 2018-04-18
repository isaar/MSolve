using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Output
{
    /// <summary>
    /// Numbers will be displayed in scientific (exponential) format. E.g -1.05E+003 
    /// </summary>
    public class ExponentialFormat: INumericFormat
    {
        /// <summary>
        /// Decimal digits to display, not counting the decimal separator (. or ,).
        /// </summary>
        public int NumDecimalDigits { get; set; } = 6;

        public string GetRealNumberFormat()
        {
            // -1.05E+003 => sign + integer + decimal separator + E + exponent sign + exponent digits + decimal digits.
            int totalWidth = NumDecimalDigits + 8;
            return "{0," + $"{totalWidth}:E{NumDecimalDigits}" + "}";
        }
    }
}
