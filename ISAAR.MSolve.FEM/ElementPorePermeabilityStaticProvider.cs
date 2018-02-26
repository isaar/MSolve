//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;
//using AnalSharp.PreProcessor.Interfaces;
//using AnalSharp.Matrices.Interfaces;
//using AnalSharp.Matrices;

//namespace AnalSharp.PreProcessor
//{
//    public class ElementPorePermeabilityStaticProvider : IElementMatrixProvider
//    {
//        private IMatrix2D<double> PorousMatrix(Element element)
//        {
//            IPorousFiniteElement elementType = (IPorousFiniteElement)element.ElementType;
//            int dofs = 0;
//            foreach (IList<DOFType> dofTypes in elementType.DOFTypes)
//                foreach (DOFType dofType in dofTypes) dofs++;
//            SymmetricMatrix2D<double> poreStiffness = new SymmetricMatrix2D<double>(dofs);
//            for (int i = 0; i < dofs; i++) poreStiffness[i, i] = 1;

//            IMatrix2D<double> permeability = elementType.PermeabilityMatrix(element);
//            int matrixRow = 0;
//            int solidRow = 0;
//            int fluidRow = 0;
//            foreach (IList<DOFType> dofTypesRow in elementType.DOFTypes)
//                foreach (DOFType dofTypeRow in dofTypesRow)
//                {
//                    int matrixCol = 0;
//                    int solidCol = 0;
//                    int fluidCol = 0;
//                    foreach (IList<DOFType> dofTypesCol in elementType.DOFTypes)
//                        foreach (DOFType dofTypeCol in dofTypesCol)
//                        {
//                            if (dofTypeCol == DOFType.Pore)
//                            {
//                                if (dofTypeRow == DOFType.Pore)
//                                    poreStiffness[matrixRow, matrixCol] = permeability[fluidRow, fluidCol];
//                                fluidCol++;
//                            }
//                            else
//                            {
//                                solidCol++;
//                            }
//                            matrixCol++;
//                        }

//                    if (dofTypeRow == DOFType.Pore)
//                        fluidRow++;
//                    else
//                        solidRow++;
//                    matrixRow++;
//                }

//            return poreStiffness;
//        }

//        #region IElementMatrixProvider Members

//        public IMatrix2D<double> Matrix(Element element)
//        {
//            if (element.ElementType is IPorousFiniteElement)
//                return PorousMatrix(element);
//            else
//            {
//                int dofs = 0;
//                foreach (DOFType[] dofTypes in element.ElementType.DOFTypes)
//                    foreach (DOFType dofType in dofTypes)
//                        dofs++;
//                SymmetricMatrix2D<double> stiffnessMatrix = new SymmetricMatrix2D<double>(dofs);
//                for (int i = 0; i < dofs; i++) stiffnessMatrix[i, i] = 1;
//                return stiffnessMatrix;
//            }
//        }

//        #endregion
//    }
//}
