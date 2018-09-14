using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.IGA.Interfaces
{
    public interface IElementMatrixProvider
    {
        IMatrix2D Matrix(Element element);
    }
}
