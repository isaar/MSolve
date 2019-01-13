using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.FEM.Providers
{
    public class ElementPoreDampingProvider_v2 : IElementMatrixProvider_v2
    {
        private readonly IElementMatrixProvider_v2 solidDampingProvider;
        private readonly double dampingCoefficient;

        public ElementPoreDampingProvider_v2(IElementMatrixProvider_v2 solidDampingProvider, double dampingCoefficient)
        {
            this.solidDampingProvider = solidDampingProvider;
            this.dampingCoefficient = dampingCoefficient;
        }

        private IMatrix PorousMatrix(IElement_v2 element)
        {
            IPorousFiniteElement_v2 elementType = (IPorousFiniteElement_v2)element.ElementType;
            int dofs = 0;
            foreach (IList<DOFType> dofTypes in elementType.DofEnumerator.GetDOFTypes(element))
                foreach (DOFType dofType in dofTypes) dofs++;
            var poreDamping = SymmetricMatrix.CreateZero(dofs);

            IMatrix damping = solidDampingProvider.Matrix(element);
            IMatrix saturation = elementType.SaturationMatrix(element);
            IMatrix coupling = elementType.CouplingMatrix(element);

            int matrixRow = 0;
            int solidRow = 0;
            int fluidRow = 0;
            foreach (IList<DOFType> dofTypesRow in elementType.DofEnumerator.GetDOFTypes(element))
                foreach (DOFType dofTypeRow in dofTypesRow)
                {
                    int matrixCol = 0;
                    int solidCol = 0;
                    int fluidCol = 0;
                    foreach (IList<DOFType> dofTypesCol in elementType.DofEnumerator.GetDOFTypes(element))
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

        public IMatrix Matrix(IElement_v2 element)
        {
            if (element.ElementType is IPorousFiniteElement_v2)
                return PorousMatrix(element);
            else
            {
                IMatrix dampingMatrix = solidDampingProvider.Matrix(element);
                dampingMatrix.Scale(dampingCoefficient);
                return dampingMatrix;

            }
        }

        #endregion
    }
}
