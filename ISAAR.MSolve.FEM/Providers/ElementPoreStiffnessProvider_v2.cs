using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.FEM.Providers
{
    public class ElementPoreStiffnessProvider_v2 : IElementMatrixProvider_v2
    {
        private readonly IElementMatrixProvider_v2 solidStiffnessProvider;
        private readonly double stiffnessCoefficient;

        public ElementPoreStiffnessProvider_v2(IElementMatrixProvider_v2 solidStiffnessProvider, double stiffnessCoefficient)
        {
            this.solidStiffnessProvider = solidStiffnessProvider;
            this.stiffnessCoefficient = stiffnessCoefficient;
        }

        private IMatrix PorousMatrix(IElement_v2 element)
        {
            IPorousFiniteElement_v2 elementType = (IPorousFiniteElement_v2)element.ElementType;
            int dofs = 0;
            foreach (IList<DOFType> dofTypes in elementType.DofEnumerator.GetDOFTypes(element))
                foreach (DOFType dofType in dofTypes) dofs++;
            var poreStiffness = SymmetricMatrix.CreateZero(dofs);

            IMatrix stiffness = solidStiffnessProvider.Matrix(element);
            IMatrix permeability = elementType.PermeabilityMatrix(element);

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

        #region IElementMatrixProvider_v2 Members

        public IMatrix Matrix(IElement_v2 element)
        {
            if (element.ElementType is IPorousFiniteElement_v2)
                return PorousMatrix(element);
            else
            {
                IMatrix stiffnessMatrix = solidStiffnessProvider.Matrix(element);
                stiffnessMatrix.Scale(stiffnessCoefficient);
                return stiffnessMatrix;
            }
        }

        #endregion
    }
}
