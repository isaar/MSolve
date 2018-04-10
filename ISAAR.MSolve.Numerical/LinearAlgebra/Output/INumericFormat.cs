using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Output
{
    /// <summary>
    /// Specifies formatting and alignment for <see cref="double"/> matrix or vector entries.
    /// </summary>
    public interface INumericFormat //TODO: Perhaps I should also specify culture
    {
        string GetRealNumberFormat();
    }
}
