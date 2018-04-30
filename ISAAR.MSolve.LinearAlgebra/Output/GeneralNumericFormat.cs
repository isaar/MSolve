using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    /// <summary>
    /// Most compact numeric format. However the matrix/vector rows will not be aligned.
    /// </summary>
    public class GeneralNumericFormat: INumericFormat
    {
        public string GetRealNumberFormat()
        {
            return "{0:G}";
        }
    }
}
