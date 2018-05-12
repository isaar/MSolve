using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using IElementMatrixProvider = ISAAR.MSolve.Discretization.Interfaces.IElementMatrixProvider;

namespace ISAAR.MSolve.FEM.Providers
{
    public class ElementPoreDampingProvider : IElementMatrixProvider
    {
        private readonly IElementMatrixProvider solidDampingProvider;
        private readonly double dampingCoefficient;

        public ElementPoreDampingProvider(IElementMatrixProvider solidDampingProvider, double dampingCoefficient)
        {
            this.solidDampingProvider = solidDampingProvider;
            this.dampingCoefficient = dampingCoefficient;
        }

        private IMatrix2D PorousMatrix(IElement element)
        {
            IPorousFiniteElement elementType = (IPorousFiniteElement)element.IElementType;
            int dofs = 0;
            foreach (IList<DOFType> dofTypes in elementType.DOFEnumerator.GetDOFTypes(element))
                foreach (DOFType dofType in dofTypes) dofs++;
            SymmetricMatrix2D poreDamping = new SymmetricMatrix2D(dofs);

            IMatrix2D damping = solidDampingProvider.Matrix(element);
            IMatrix2D saturation = elementType.SaturationMatrix(element);
            IMatrix2D coupling = elementType.CouplingMatrix(element);

            int matrixRow = 0;
            int solidRow = 0;
            int fluidRow = 0;
            foreach (IList<DOFType> dofTypesRow in elementType.DOFEnumerator.GetDOFTypes(element))
                foreach (DOFType dofTypeRow in dofTypesRow)
                {
                    int matrixCol = 0;
                    int solidCol = 0;
                    int fluidCol = 0;
                    foreach (IList<DOFType> dofTypesCol in elementType.DOFEnumerator.GetDOFTypes(element))
                        foreach (DOFType dofTypeCol in dofTypesCol)
                        {
                            if (dofTypeCol == DOFType.Pore)
                            {
                                if (dofTypeRow == DOFType.Pore)
                                    poreDamping[matrixRow, matrixCol] = -saturation[fluidRow, fluidCol];
                                else
                                    poreDamping[matrixRow, matrixCol] = coupling[fluidCol, solidRow];
                                fluidCol++;
                            }
                            else
                            {
                                if (dofTypeRow != DOFType.Pore)
                                    poreDamping[matrixRow, matrixCol] = damping[solidRow, solidCol] * dampingCoefficient;
                                else
                                    poreDamping[matrixRow, matrixCol] = coupling[fluidRow, solidCol];
                                solidCol++;
                            }
                            matrixCol++;
                        }

                    if (dofTypeRow == DOFType.Pore)
                        fluidRow++;
                    else
                        solidRow++;
                    matrixRow++;
                }

            return poreDamping;
        }

        #region IElementMatrixProvider Members

        public IMatrix2D Matrix(IElement element)
        {
            if (element.IElementType is IPorousFiniteElement)
                return PorousMatrix(element);
            else
            {
                IMatrix2D dampingMatrix = solidDampingProvider.Matrix(element);
                dampingMatrix.Scale(dampingCoefficient);
                return dampingMatrix;

            }
        }

        #endregion
    }
}
