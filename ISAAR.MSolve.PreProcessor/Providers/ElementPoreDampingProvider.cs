using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.PreProcessor.Providers
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

        private IMatrix2D<double> PorousMatrix(Element element)
        {
            IPorousFiniteElement elementType = (IPorousFiniteElement)element.ElementType;
            int dofs = 0;
            foreach (IList<DOFType> dofTypes in elementType.DOFEnumerator.GetDOFTypes(element))
                foreach (DOFType dofType in dofTypes) dofs++;
            SymmetricMatrix2D<double> poreDamping = new SymmetricMatrix2D<double>(dofs);

            IMatrix2D<double> damping = solidDampingProvider.Matrix(element);
            IMatrix2D<double> saturation = elementType.SaturationMatrix(element);
            IMatrix2D<double> coupling = elementType.CouplingMatrix(element);

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

        public IMatrix2D<double> Matrix(Element element)
        {
            if (element.ElementType is IPorousFiniteElement)
                return PorousMatrix(element);
            else
            {
                IMatrix2D<double> dampingMatrix = solidDampingProvider.Matrix(element);
                dampingMatrix.Scale(dampingCoefficient);
                return dampingMatrix;

            }
        }

        #endregion
    }
}
