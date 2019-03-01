using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.MultiscaleAnalysisMerge
{
    /// <summary>
    /// Element stiffness matrix provider for simultaneous calculation of global stiffness matrix macroscopic variables in multiscale FE2 scheme
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class StiffnessProviderSimu_v2 : IElementMatrixProvider_v2
    {
        #region IElementMatrixProvider Members
        private SubdomainCalculationsAndAssembly host;

        public StiffnessProviderSimu_v2(SubdomainCalculationsAndAssembly host)
        {
            this.host = host;
        }

        public IMatrix Matrix(IElement_v2 element) //TODOGer IMatrix2D will be changed to Matrix etc.
        {
            var elementMatrix = element.ElementType.StiffnessMatrix(element);
            host.UpdateVectors_v2(element, elementMatrix);
            return elementMatrix;
        }

        #endregion
    }
}






    
        //#region IElementMatrixProvider Members

        //public IMatrix2D Matrix(IElement element)
        //{
        //    return element.IElementType.StiffnessMatrix(element);
        //}

        //#endregion
    
