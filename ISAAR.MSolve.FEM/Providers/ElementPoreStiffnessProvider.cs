using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using IElementMatrixProvider = ISAAR.MSolve.Discretization.Interfaces.IElementMatrixProvider;

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

        private IMatrix2D PorousMatrix(IElement element)
        {
            IPorousFiniteElement elementType = (IPorousFiniteElement)element.IElementType;
            int dofs = 0;
            foreach (IList<DOFType> dofTypes in elementType.DOFEnumerator.GetDOFTypes(element))
                foreach (DOFType dofType in dofTypes) dofs++;
            SymmetricMatrix2D poreStiffness = new SymmetricMatrix2D(dofs);

            IMatrix2D stiffness = solidStiffnessProvider.Matrix(element);
            IMatrix2D permeability = elementType.PermeabilityMatrix(element);

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
                                    // H correction
                                    poreStiffness[matrixRow, matrixCol] = -permeability[fluidRow, fluidCol];
                                    //poreStiffness[matrixRow, matrixCol] = permeability[fluidRow, fluidCol];
                                fluidCol++;
                            }
                            else
                            {
                                if (dofTypeRow != DOFType.Pore)
                                    poreStiffness[matrixRow, matrixCol] = stiffness[solidRow, solidCol] * stiffnessCoefficient;
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

            return poreStiffness;
        }

        #region IElementMatrixProvider Members

        public IMatrix2D Matrix(IElement element)
        {
            if (element.IElementType is IPorousFiniteElement)
                return PorousMatrix(element);
            else
            {
                IMatrix2D stiffnessMatrix = solidStiffnessProvider.Matrix(element);
                stiffnessMatrix.Scale(stiffnessCoefficient);
                return stiffnessMatrix;
            }
        }

        #endregion
    }
}
