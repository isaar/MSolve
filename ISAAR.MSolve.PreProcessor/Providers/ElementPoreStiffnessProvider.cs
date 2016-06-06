using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.PreProcessor.Providers
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

        private IMatrix2D<double> PorousMatrix(Element element)
        {
            IPorousFiniteElement elementType = (IPorousFiniteElement)element.ElementType;
            int dofs = 0;
            foreach (IList<DOFType> dofTypes in elementType.DOFEnumerator.GetDOFTypes(element))
                foreach (DOFType dofType in dofTypes) dofs++;
            SymmetricMatrix2D<double> poreStiffness = new SymmetricMatrix2D<double>(dofs);

            IMatrix2D<double> stiffness = solidStiffnessProvider.Matrix(element);
            IMatrix2D<double> permeability = elementType.PermeabilityMatrix(element);

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

        public IMatrix2D<double> Matrix(Element element)
        {
            if (element.ElementType is IPorousFiniteElement)
                return PorousMatrix(element);
            else
            {
                IMatrix2D<double> stiffnessMatrix = solidStiffnessProvider.Matrix(element);
                stiffnessMatrix.Scale(stiffnessCoefficient);
                return stiffnessMatrix;
            }
        }

        #endregion
    }
}
