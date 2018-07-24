namespace ISAAR.MSolve.LinearAlgebra.Output.Formatting
{
    /// <summary>
    /// Most compact numeric format. Output numbers will be displayed in the default precision. However the matrix rows and
    /// vertical vector entries will not be aligned.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class GeneralNumericFormat: INumericFormat
    {
        /// <summary>
        /// Returns a tring that will be used as a template for formatting numbers during output.
        /// </summary>
        public string GetRealNumberFormat()
        {
            return "{0:G}";
        }
    }
}
