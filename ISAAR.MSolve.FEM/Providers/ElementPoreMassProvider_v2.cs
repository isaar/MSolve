using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.FEM.Providers
{
    public class ElementPoreMassProvider_v2 : IElementMatrixProvider_v2
    {
        private readonly IElementMatrixProvider_v2 solidMassProvider;
        private readonly double massCoefficient;

        public ElementPoreMassProvider_v2(IElementMatrixProvider_v2 solidMassProvider, double massCoefficient)
        {
            this.solidMassProvider = solidMassProvider;
            this.massCoefficient = massCoefficient;
        }

        private IMatrix PorousMatrix(IElement_v2 element)
        {
            IPorousFiniteElement_v2 elementType = (IPorousFiniteElement_v2)element.ElementType;
            int dofs = 0;
            foreach (IList<DOFType> dofTypes in elementType.DofEnumerator.GetDOFTypes(element))
                foreach (DOFType dofType in dofTypes) dofs++;
            var poreMass = SymmetricMatrix.CreateZero(dofs);

            IMatrix mass = solidMassProvider.Matrix(element);

            int matrixRow = 0;
            int solidRow = 0;
            foreach (IList<DOFType> dofTypesRow in elementType.DofEnumerator.GetDOFTypes(element))
                foreach (DOFType dofTypeRow in dofTypesRow)
                {
                    int matrixCol = 0;
                    int solidCol = 0;
                    foreach (IList<DOFType> dofTypesCol in elementType.DofEnumerator.GetDOFTypes(element))
                        foreach (DOFType dofTypeCol in dofTypesCol)
                        {
                            if (dofTypeCol == DOFType.Pore)
                            {
                            }
                            else
                            {
                                if (dofTypeRow != DOFType.Pore)
                                    poreMass[matrixRow, matrixCol] = mass[solidRow, solidCol] * massCoefficient;
                                solidCol++;
                            }
                            matrixCol++;
                        }

                    if (dofTypeRow != DOFType.Pore) solidRow++;
                    matrixRow++;
                }

            return poreMass;
        }

        #region IElementMatrixProvider Members

        public IMatrix Matrix(IElement_v2 element)
        {
            if (element.ElementType is IPorousFiniteElement_v2)
                return PorousMatrix(element);
            else
            {
                IMatrix massMatrix = solidMassProvider.Matrix(element);
                massMatrix.Scale(massCoefficient);
                return massMatrix;
            }
        }

        #endregion
    }
}
