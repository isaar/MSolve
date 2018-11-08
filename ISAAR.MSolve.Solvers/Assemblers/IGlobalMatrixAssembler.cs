//using System;
//using System.Collections.Generic;
//using System.Text;
//using ISAAR.MSolve.Discretization.Interfaces;
//using ISAAR.MSolve.LinearAlgebra.Matrices;
//using ISAAR.MSolve.Solvers.Ordering;

////TODO: not sure this interface is required
//namespace ISAAR.MSolve.Solvers.Assemblers
//{
//    public interface IGlobalMatrixAssembler<TMatrix>
//        where TMatrix: IMatrix
//    {
//        TMatrix BuildGlobalMatrix(IEnumerable<IElement> elements, FreeDofOrderer dofOrderer, 
//            IElementMatrixProvider matrixProvider);

//        TMatrix BuildGlobalMatrix(ISubdomain subdomain, IElementMatrixProvider matrixProvider);
//    }
//}
