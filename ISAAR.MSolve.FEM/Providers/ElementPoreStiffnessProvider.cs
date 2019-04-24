using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.FEM.Providers
{
    public class ElementPoreStiffnessProvider : IElementMatrixProvider
    {
        private readonly IElementMatrixProvider solidStiffnessProvider;
        private readonly double stiffnessCoefficient;

        public ElementPoreStiffnessProvider(IElementMatrixProvider solidStiffnessProvider, double stiffnessCoefficient)
        {
            this.solidStiffnessProvider = solidStiffnessProvider;
            this.stiffnessCoefficient = stiffnessCoefficient;
        }

        private IMatrix PorousMatrix(IElement element)
        {
            IPorousFiniteElement elementType = (IPorousFiniteElement)element.ElementType;
            int dofs = 0;
            foreach (IList<IDofType> dofTypes in elementType.DofEnumerator.GetDofTypesForMatrixAssembly(element))
                foreach (IDofType dofType in dofTypes) dofs++;
            var poreStiffness = SymmetricMatrix.CreateZero(dofs);

            IMatrix stiffness = solidStiffnessProvider.Matrix(element);
            IMatrix permeability = elementType.PermeabilityMatrix(element);

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
                                    // H correction
                                    poreStiffness[matrixRow, matrixCol] = -permeability[fluidRow, fluidCol];
                                    //poreStiffness[matrixRow, matrixCol] = permeability[fluidRow, fluidCol];
                                fluidCol++;
                            }
                            else
                            {
                                if (dofTypeRow != PorousMediaDof.Pressure)
                                    poreStiffness[matrixRow, matrixCol] = stiffness[solidRow, solidCol] * stiffnessCoefficient;
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

            return poreStiffness;
        }

        #region IElementMatrixProvider Members

        public IMatrix Matrix(IElement element)
        {
            if (element.ElementType is IPorousFiniteElement)
                return PorousMatrix(element);
            else
            {
                IMatrix stiffnessMatrix = solidStiffnessProvider.Matrix(element);
                stiffnessMatrix.ScaleIntoThis(stiffnessCoefficient);
                return stiffnessMatrix;
            }
        }

        #endregion
    }
}
