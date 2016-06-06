using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.PreProcessor.Interfaces
{
    public interface IElementMatrixProvider
    {
        IMatrix2D<double> Matrix(Element element);
    }
}
