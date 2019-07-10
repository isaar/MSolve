using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

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

        private IMatrix PorousMatrix(IElement element)
        {
            IPorousFiniteElement elementType = (IPorousFiniteElement)element.ElementType;
            int dofs = 0;
            foreach (IList<IDofType> dofTypes in elementType.DofEnumerator.GetDofTypesForMatrixAssembly(element))
                foreach (IDofType dofType in dofTypes) dofs++;
            var poreDamping = SymmetricMatrix.CreateZero(dofs);

            IMatrix damping = solidDampingProvider.Matrix(element);
            IMatrix saturation = elementType.SaturationMatrix(element);
            IMatrix coupling = elementType.CouplingMatrix(element);

            int matrixRow = 0;
            int solidRow = 0;
            int fluidRow = 0;
            foreach (IList<IDofType> dofTypesRow in elementType.DofEnumerator.GetDofTypesForMatrixAssembly(element))
                foreach (IDofType dofTypeRow in dofTypesRow)
                {
                    int matrixCol = 0;
                    int solidCol = 0;
                    int fluidCol = 0;
                    foreach (IList<IDofType> dofTypesCol in elementType.DofEnumerator.GetDofTypesForMatrixAssembly(element))
                        foreach (IDofType dofTypeCol in dofTypesCol)
                        {
                            if (dofTypeCol == PorousMediaDof.Pressure)
                            {
                                if (dofTypeRow == PorousMediaDof.Pressure)
                                    poreDamping[matrixRow, matrixCol] = -saturation[fluidRow, fluidCol];
                                else
                                    poreDamping[matrixRow, matrixCol] = coupling[fluidCol, solidRow];
                                fluidCol++;
                            }
                            else
                            {
                                if (dofTypeRow != PorousMediaDof.Pressure)
                                    poreDamping[matrixRow, matrixCol] = damping[solidRow, solidCol] * dampingCoefficient;
                                else
                                    poreDamping[matrixRow, matrixCol] = coupling[fluidRow, solidCol];
                                solidCol++;
                            }
                            matrixCol++;
                        }

                    if (dofTypeRow == PorousMediaDof.Pressure)
                        fluidRow++;
                    else
                        solidRow++;
                    matrixRow++;
                }

            return poreDamping;
        }

        #region IElementMatrixProvider Members

        public IMatrix Matrix(IElement element)
        {
            if (element.ElementType is IPorousFiniteElement)
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
