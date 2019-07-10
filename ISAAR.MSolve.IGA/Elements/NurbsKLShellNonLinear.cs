using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA.Elements
{
    public class NurbsKLShellNonLinear:Element,IStructuralIsogeometricElement
    {
        protected readonly static IDofType[] controlPointDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
        protected IDofType[][] dofTypes;
        protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
        public CellType CellType { get; } = CellType.Unknown;

        public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;
        public IElementDofEnumerator DofEnumerator
        {
            get => dofEnumerator;
            set => this.dofEnumerator = value;
        }
        public bool MaterialModified { get; }
        public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, NeumannBoundaryCondition neumann)
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann)
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure)
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure)
        {
            throw new NotImplementedException();
        }

        public void ResetMaterialModified()
        {
            throw new NotImplementedException();
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[,] CalculateDisplacementsForPostProcessing(Element element, double[,] localDisplacements)
        {
            var nurbsElement = (NURBSKirchhoffLoveShellElement)element;
            var knotParametricCoordinatesKsi = Vector.CreateFromArray(new double[] { element.Knots[0].Ksi, element.Knots[2].Ksi });
            var knotParametricCoordinatesHeta = Vector.CreateFromArray(new double[] { element.Knots[0].Heta, element.Knots[1].Heta });
            NURBS2D nurbs = new NURBS2D(nurbsElement, nurbsElement.ControlPoints, knotParametricCoordinatesKsi, knotParametricCoordinatesHeta);
            var knotDisplacements = new double[4, 3];
            var paraviewKnotRenumbering = new int[] { 0, 3, 1, 2 };
            for (int j = 0; j < element.Knots.Count; j++)
            {
                for (int i = 0; i < element.ControlPoints.Count; i++)
                {
                    knotDisplacements[paraviewKnotRenumbering[j], 0] += nurbs.NurbsValues[i, j] * localDisplacements[i, 0];
                    knotDisplacements[paraviewKnotRenumbering[j], 1] += nurbs.NurbsValues[i, j] * localDisplacements[i, 1];
                    knotDisplacements[paraviewKnotRenumbering[j], 2] += nurbs.NurbsValues[i, j] * localDisplacements[i, 2];
                }
            }

            return knotDisplacements;
        }

        public void ClearMaterialState()
        {
            throw new NotImplementedException();
        }

        private Vector[] initialSurfaceBasisVectors1;
        private Vector[] initialSurfaceBasisVectors2;
        private Vector[] initialSurfaceBasisVectors3;

        private Vector[] initialSurfaceBasisVectorDerivative1;
        private Vector[] initialSurfaceBasisVectorDerivative2;
        private Vector[] initialSurfaceBasisVectorDerivative12;

        private bool isInitialized;
        public IMatrix StiffnessMatrix(IElement element)
        {
            var shellElement = (NurbsKLShellNonLinear)element;

            IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(shellElement);
            Matrix stiffnessMatrixElement = Matrix.CreateZero(shellElement.ControlPointsDictionary.Count * 3, shellElement.ControlPointsDictionary.Count * 3);

            NURBS2D nurbs = new NURBS2D(shellElement, shellElement.ControlPoints);

            if (!isInitialized)
            {
                CalculateInitialConfigurationData(shellElement, nurbs,gaussPoints);
                //var localTotalDisplacements = new double[shellElement.ControlPoints.Count*3];
                //var aaaaa = UpdateCoordinateData(localTotalDisplacements);
                //CalculateStrains(localTotalDisplacements, element, aaaaa);
                isInitialized = true;
            }

            for (int j = 0; j < gaussPoints.Count; j++)
            {
                var jacobianMatrix = CalculateJacobian(shellElement, nurbs, j);

                var hessianMatrix = CalculateHessian(shellElement, nurbs, j);
                var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

                var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

                var surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);

                var J1 = surfaceBasisVector3.Norm2();
                surfaceBasisVector3.ScaleIntoThis(1 / J1);

                var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
                var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
                var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

                Matrix constitutiveMatrix = CalculateConstitutiveMatrix(shellElement, surfaceBasisVector1, surfaceBasisVector2);

                var Bmembrane = CalculateMembraneDeformationMatrix(nurbs, j, surfaceBasisVector1, surfaceBasisVector2, shellElement);

                var Bbending = CalculateBendingDeformationMatrix(surfaceBasisVector3, nurbs, j, surfaceBasisVector2, surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2, surfaceBasisVectorDerivative12, shellElement);

                double membraneStiffness = ((IIsotropicContinuumMaterial2D)shellElement.Patch.Material).YoungModulus * shellElement.Patch.Thickness /
                                           (1 - Math.Pow(((IIsotropicContinuumMaterial2D)shellElement.Patch.Material).PoissonRatio, 2));

                var Kmembrane = Bmembrane.Transpose() * constitutiveMatrix * Bmembrane * membraneStiffness * J1 *
                                gaussPoints[j].WeightFactor;

                var KmembraneNL = CalculateKmembraneNL(shellElement, membraneStiffness, constitutiveMatrix, surfaceBasisVector1, j, surfaceBasisVector2, nurbs);

                

                double bendingStiffness = ((IIsotropicContinuumMaterial2D)shellElement.Patch.Material).YoungModulus * Math.Pow(shellElement.Patch.Thickness, 3) /
                                          12 / (1 - Math.Pow(((IIsotropicContinuumMaterial2D)shellElement.Patch.Material).PoissonRatio, 2));


                var Kbending = Bbending.Transpose() * constitutiveMatrix * Bbending * bendingStiffness * J1 *
                               gaussPoints[j].WeightFactor;


                var KbendingNL = Matrix.CreateZero(shellElement.ControlPoints.Count * 3, shellElement.ControlPoints.Count * 3);
                var bendingStrains = bendingStiffness * constitutiveMatrix * Vector.CreateFromArray(new double[]
                {
                    0.5*(surfaceBasisVectorDerivative1*surfaceBasisVectorDerivative1-initialSurfaceBasisVectorDerivative1[j]*initialSurfaceBasisVectorDerivative1[j]),
                    0.5*(surfaceBasisVectorDerivative2*surfaceBasisVectorDerivative2-initialSurfaceBasisVectorDerivative2[j]*initialSurfaceBasisVectorDerivative2[j]),
                    surfaceBasisVectorDerivative1*surfaceBasisVectorDerivative2-initialSurfaceBasisVectorDerivative1[j]*initialSurfaceBasisVectorDerivative2[j]
                });

                for (int i = 0; i < shellElement.ControlPoints.Count; i++)
                {
                    for (int k = 0; k < shellElement.ControlPoints.Count; k++)
                    {
                        var a11r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsi[i, j]);
                        var a22r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueHeta[i, j]);
                        var a12r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsiHeta[i, j]);


                        var a11s = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsi[k, j]);
                        var a22s = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueHeta[k, j]);
                        var a12s = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsiHeta[k, j]);

                        var a3r = CalculateA3r(nurbs, i, j, surfaceBasisVector2, surfaceBasisVector1);
                        var a3s = CalculateA3r(nurbs, k, j, surfaceBasisVector2, surfaceBasisVector1);

                        var a1r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesKsi[i, j]);
                        var a2s = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesHeta[k, j]);

                        //var a1
                    }
                }

                stiffnessMatrixElement.AddIntoThis(Kmembrane);
                stiffnessMatrixElement.AddIntoThis(KmembraneNL);

                stiffnessMatrixElement.AddIntoThis(Kbending);
            }
            return stiffnessMatrixElement;
        }

        private Matrix3by3 CalculateA3r(NURBS2D nurbs, int i, int j, Vector surfaceBasisVector2, Vector surfaceBasisVector1)
        {
            var aux1 = Vector.CreateFromArray(new double[] {nurbs.NurbsDerivativeValuesKsi[i, j], 0, 0})
                .CrossProduct(surfaceBasisVector2);
            var aux2 = Vector.CreateFromArray(new double[] {0, nurbs.NurbsDerivativeValuesKsi[i, j], 0})
                .CrossProduct(surfaceBasisVector2);
            var aux3 = Vector.CreateFromArray(new double[] {0, 0, nurbs.NurbsDerivativeValuesKsi[i, j]})
                .CrossProduct(surfaceBasisVector2);

            var aux4 = surfaceBasisVector1.CrossProduct(Vector.CreateFromArray(new double[]
                {nurbs.NurbsDerivativeValuesHeta[i, j], 0, 0}));
            var aux5 = surfaceBasisVector1.CrossProduct(Vector.CreateFromArray(new double[]
                {0, nurbs.NurbsDerivativeValuesHeta[i, j], 0}));
            var aux6 = surfaceBasisVector1.CrossProduct(Vector.CreateFromArray(new double[]
                {0, 0, nurbs.NurbsDerivativeValuesHeta[i, j]}));

            var a3r = Matrix3by3.CreateZero();
            a3r[0, 0] = aux1[0] + aux4[0];
            a3r[0, 0] = aux1[1] + aux4[1];
            a3r[0, 0] = aux1[2] + aux4[2];

            a3r[0, 0] = aux2[0] + aux5[0];
            a3r[0, 0] = aux2[1] + aux5[1];
            a3r[0, 0] = aux2[2] + aux5[2];

            a3r[0, 0] = aux3[0] + aux6[0];
            a3r[0, 0] = aux3[1] + aux6[1];
            a3r[0, 0] = aux3[2] + aux6[2];
            return a3r;
        }

        private Matrix CalculateKmembraneNL(NurbsKLShellNonLinear shellElement, double membraneStiffness,
            Matrix constitutiveMatrix, Vector surfaceBasisVector1, int j, Vector surfaceBasisVector2, NURBS2D nurbs)
        {
            var KmembraneNL = Matrix.CreateZero(shellElement.ControlPoints.Count * 3, shellElement.ControlPoints.Count * 3);
            var membraneStrains = membraneStiffness * constitutiveMatrix * Vector.CreateFromArray(new double[]
            {
                0.5 * (surfaceBasisVector1 * surfaceBasisVector1 -
                       initialSurfaceBasisVectors1[j] * initialSurfaceBasisVectors1[j]),
                0.5 * (surfaceBasisVector2 * surfaceBasisVector2 -
                       initialSurfaceBasisVectors2[j] * initialSurfaceBasisVectors2[j]),
                surfaceBasisVector1 * surfaceBasisVector2 -
                initialSurfaceBasisVectors1[j] * initialSurfaceBasisVectors2[j]
            });

            for (int i = 0; i < shellElement.ControlPoints.Count; i++)
            {
                for (int k = 0; k < shellElement.ControlPoints.Count; k++)
                {
                    var a1r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesKsi[i, j]);
                    var a1s = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesKsi[k, j]);
                    var a2r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesHeta[i, j]);
                    var a2s = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesHeta[k, j]);

                    var klocal = membraneStrains[0] * a1r * a1s + membraneStrains[1] * a2r * a2s +
                                 membraneStrains[2] * (a1r * a2s + a1s * a2r);

                    for (int l = i * 3; l < i * 3 + 3; l++)
                    for (int m = k * 3; m < k * 3 + 3; m++)
                        KmembraneNL[l, m] += klocal[l / 3, m / 3];
                }
            }

            return KmembraneNL;
        }

        private void CalculateInitialConfigurationData(NurbsKLShellNonLinear shellElement, NURBS2D nurbs, IList<GaussLegendrePoint3D> gaussPoints)
        {
            var numberOfGP = gaussPoints.Count;
            initialSurfaceBasisVectors1= new Vector[numberOfGP];
            initialSurfaceBasisVectors2 = new Vector[numberOfGP];
            initialSurfaceBasisVectors3 = new Vector[numberOfGP];
            initialSurfaceBasisVectorDerivative1 = new Vector[numberOfGP];
            initialSurfaceBasisVectorDerivative2 = new Vector[numberOfGP];
            initialSurfaceBasisVectorDerivative12 = new Vector[numberOfGP];

            for (int j = 0; j < gaussPoints.Count; j++)
            {
                var jacobianMatrix = CalculateJacobian(shellElement, nurbs, j);

                var hessianMatrix = CalculateHessian(shellElement, nurbs, j);
                initialSurfaceBasisVectors1[j] = CalculateSurfaceBasisVector1(jacobianMatrix, 0);
                initialSurfaceBasisVectors2[j] = CalculateSurfaceBasisVector1(jacobianMatrix, 1);
                initialSurfaceBasisVectors3[j] = initialSurfaceBasisVectors1[j].CrossProduct(initialSurfaceBasisVectors2[j]);

                initialSurfaceBasisVectorDerivative1[j] = CalculateSurfaceBasisVector1(hessianMatrix, 0);
                initialSurfaceBasisVectorDerivative2[j] = CalculateSurfaceBasisVector1(hessianMatrix, 1);
                initialSurfaceBasisVectorDerivative12[j] = CalculateSurfaceBasisVector1(hessianMatrix, 2);
            }
        }

        private Matrix CalculateConstitutiveMatrix(NurbsKLShellNonLinear element, Vector surfaceBasisVector1, Vector surfaceBasisVector2)
        {
            var auxMatrix1 = Matrix.CreateZero(2, 2);
            auxMatrix1[0, 0] = surfaceBasisVector1.DotProduct(surfaceBasisVector1);
            auxMatrix1[0, 1] = surfaceBasisVector1.DotProduct(surfaceBasisVector2);
            auxMatrix1[1, 0] = surfaceBasisVector2.DotProduct(surfaceBasisVector1);
            auxMatrix1[1, 1] = surfaceBasisVector2.DotProduct(surfaceBasisVector2);
            (Matrix inverse, double det) = auxMatrix1.InvertAndDeterminant();

            var material = ((IContinuumMaterial2D)element.Patch.Material);
            var constitutiveMatrix = Matrix.CreateFromArray(new double[3, 3]
            {
                {
                    inverse[0,0]*inverse[0,0],
                    material.PoissonRatio*inverse[0,0]*inverse[1,1]+(1-material.PoissonRatio)*inverse[1,0]*inverse[1,0],
                    inverse[0,0]*inverse[1,0]
                },
                {
                    material.PoissonRatio*inverse[0,0]*inverse[1,1]+(1-material.PoissonRatio)*inverse[1,0]*inverse[1,0],
                    inverse[1,1]*inverse[1,1],
                    inverse[1,1]*inverse[1,0]
                },
                {
                    inverse[0,0]*inverse[1,0],
                    inverse[1,1]*inverse[1,0],
                    0.5*(1-material.PoissonRatio)*inverse[0,0]*inverse[1,1]+(1+material.PoissonRatio)*inverse[1,0]*inverse[1,0]
                },
            });
            return constitutiveMatrix;
        }

        private Matrix CalculateBendingDeformationMatrix(Vector surfaceBasisVector3, NURBS2D nurbs, int j,
            Vector surfaceBasisVector2, Vector surfaceBasisVectorDerivative1, Vector surfaceBasisVector1, double J1,
            Vector surfaceBasisVectorDerivative2, Vector surfaceBasisVectorDerivative12, NurbsKLShellNonLinear element)
        {
            Matrix Bbending = Matrix.CreateZero(3, element.ControlPoints.Count * 3);
            for (int column = 0; column < element.ControlPoints.Count * 3; column += 3)
            {
                #region BI1

                var BI1 = surfaceBasisVector3.CrossProduct(surfaceBasisVector3);
                BI1.ScaleIntoThis(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
                var auxVector = surfaceBasisVector2.CrossProduct(surfaceBasisVector3);
                auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
                BI1.AddIntoThis(auxVector);
                BI1.ScaleIntoThis(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative1));
                auxVector = surfaceBasisVector1.CrossProduct(surfaceBasisVectorDerivative1);
                auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
                BI1.AddIntoThis(auxVector);
                BI1.ScaleIntoThis(1 / J1);
                auxVector[0] = surfaceBasisVector3[0];
                auxVector[1] = surfaceBasisVector3[1];
                auxVector[2] = surfaceBasisVector3[2];
                auxVector.ScaleIntoThis(-nurbs.NurbsSecondDerivativeValueKsi[column / 3, j]);
                BI1.AddIntoThis(auxVector);

                #endregion

                #region BI2

                Vector BI2 = surfaceBasisVector3.CrossProduct(surfaceBasisVector3);
                BI2.ScaleIntoThis(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
                auxVector = surfaceBasisVector2.CrossProduct(surfaceBasisVector3);
                auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
                BI2.AddIntoThis(auxVector);
                BI2.ScaleIntoThis(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative2));
                auxVector = surfaceBasisVector1.CrossProduct(surfaceBasisVectorDerivative2);
                auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
                BI2.AddIntoThis(auxVector);
                auxVector = surfaceBasisVectorDerivative2.CrossProduct(surfaceBasisVector2);
                auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
                BI2.AddIntoThis(auxVector);
                BI2.ScaleIntoThis(1 / J1);
                auxVector[0] = surfaceBasisVector3[0];
                auxVector[1] = surfaceBasisVector3[1];
                auxVector[2] = surfaceBasisVector3[2];
                auxVector.ScaleIntoThis(-nurbs.NurbsSecondDerivativeValueHeta[column / 3, j]);
                BI2.AddIntoThis(auxVector);

                #endregion

                #region BI3

                Vector BI3 = surfaceBasisVector3.CrossProduct(surfaceBasisVector3);
                BI3.ScaleIntoThis(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
                auxVector = surfaceBasisVector2.CrossProduct(surfaceBasisVector3);
                auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
                BI3.AddIntoThis(auxVector);
                BI3.ScaleIntoThis(surfaceBasisVector3.DotProduct(surfaceBasisVectorDerivative12));
                auxVector = surfaceBasisVector1.CrossProduct(surfaceBasisVectorDerivative12);
                auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesHeta[column / 3, j]);
                BI3.AddIntoThis(auxVector);
                auxVector = surfaceBasisVectorDerivative2.CrossProduct(surfaceBasisVector2);
                auxVector.ScaleIntoThis(nurbs.NurbsDerivativeValuesKsi[column / 3, j]);
                BI3.AddIntoThis(auxVector);
                BI3.ScaleIntoThis(1 / J1);
                auxVector[0] = surfaceBasisVector3[0];
                auxVector[1] = surfaceBasisVector3[1];
                auxVector[2] = surfaceBasisVector3[2];
                auxVector.ScaleIntoThis(-nurbs.NurbsSecondDerivativeValueKsiHeta[column / 3, j]);
                BI3.AddIntoThis(auxVector);

                #endregion

                Bbending[0, column] = BI1[0];
                Bbending[0, column + 1] = BI1[1];
                Bbending[0, column + 2] = BI1[2];

                Bbending[1, column] = BI2[0];
                Bbending[1, column + 1] = BI2[1];
                Bbending[1, column + 2] = BI2[2];

                Bbending[2, column] = 2 * BI3[0];
                Bbending[2, column + 1] = 2 * BI3[1];
                Bbending[2, column + 2] = 2 * BI3[2];
            }

            return Bbending;
        }

        private Matrix CalculateMembraneDeformationMatrix(NURBS2D nurbs, int j, Vector surfaceBasisVector1,
            Vector surfaceBasisVector2, NurbsKLShellNonLinear element)
        {
            Matrix dRIa = Matrix.CreateZero(3, element.ControlPoints.Count);
            for (int i = 0; i < element.ControlPoints.Count; i++)
            {
                for (int m = 0; m < 3; m++)
                {
                    dRIa[m, i] = nurbs.NurbsDerivativeValuesHeta[i, j] * surfaceBasisVector1[m] +
                                 nurbs.NurbsDerivativeValuesKsi[i, j] * surfaceBasisVector2[m];
                }
            }

            Matrix Bmembrane = Matrix.CreateZero(3, element.ControlPoints.Count * 3);
            for (int column = 0; column < element.ControlPoints.Count * 3; column += 3)
            {
                Bmembrane[0, column] = nurbs.NurbsDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[0];
                Bmembrane[0, column + 1] = nurbs.NurbsDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[1];
                Bmembrane[0, column + 2] = nurbs.NurbsDerivativeValuesKsi[column / 3, j] * surfaceBasisVector1[2];

                Bmembrane[1, column] = nurbs.NurbsDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[0];
                Bmembrane[1, column + 1] = nurbs.NurbsDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[1];
                Bmembrane[1, column + 2] = nurbs.NurbsDerivativeValuesHeta[column / 3, j] * surfaceBasisVector2[2];

                Bmembrane[2, column] = dRIa[0, column / 3];
                Bmembrane[2, column + 1] = dRIa[1, column / 3];
                Bmembrane[2, column + 2] = dRIa[2, column / 3];
            }

            return Bmembrane;
        }

        private static Vector CalculateSurfaceBasisVector1(Matrix Matrix, int row)
        {
            Vector surfaceBasisVector1 = Vector.CreateZero(3);
            surfaceBasisVector1[0] = Matrix[row, 0];
            surfaceBasisVector1[1] = Matrix[row, 1];
            surfaceBasisVector1[2] = Matrix[row, 2];
            return surfaceBasisVector1;
        }

        private static Matrix CalculateHessian(NurbsKLShellNonLinear shellElement, NURBS2D nurbs, int j)
        {
            Matrix hessianMatrix = Matrix.CreateZero(3, 3);
            for (int k = 0; k < shellElement.ControlPoints.Count; k++)
            {
                hessianMatrix[0, 0] += nurbs.NurbsSecondDerivativeValueKsi[k, j] * shellElement.ControlPoints[k].X;
                hessianMatrix[0, 1] += nurbs.NurbsSecondDerivativeValueKsi[k, j] * shellElement.ControlPoints[k].Y;
                hessianMatrix[0, 2] += nurbs.NurbsSecondDerivativeValueKsi[k, j] * shellElement.ControlPoints[k].Z;
                hessianMatrix[1, 0] += nurbs.NurbsSecondDerivativeValueHeta[k, j] * shellElement.ControlPoints[k].X;
                hessianMatrix[1, 1] += nurbs.NurbsSecondDerivativeValueHeta[k, j] * shellElement.ControlPoints[k].Y;
                hessianMatrix[1, 2] += nurbs.NurbsSecondDerivativeValueHeta[k, j] * shellElement.ControlPoints[k].Z;
                hessianMatrix[2, 0] += nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * shellElement.ControlPoints[k].X;
                hessianMatrix[2, 1] += nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * shellElement.ControlPoints[k].Y;
                hessianMatrix[2, 2] += nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * shellElement.ControlPoints[k].Z;
            }

            return hessianMatrix;
        }

        private static Matrix CalculateJacobian(NurbsKLShellNonLinear shellElement, NURBS2D nurbs, int j)
        {
            Matrix jacobianMatrix = Matrix.CreateZero(2, 3);
            for (int k = 0; k < shellElement.ControlPoints.Count; k++)
            {
                jacobianMatrix[0, 0] += nurbs.NurbsDerivativeValuesKsi[k, j] * shellElement.ControlPoints[k].X;
                jacobianMatrix[0, 1] += nurbs.NurbsDerivativeValuesKsi[k, j] * shellElement.ControlPoints[k].Y;
                jacobianMatrix[0, 2] += nurbs.NurbsDerivativeValuesKsi[k, j] * shellElement.ControlPoints[k].Z;
                jacobianMatrix[1, 0] += nurbs.NurbsDerivativeValuesHeta[k, j] * shellElement.ControlPoints[k].X;
                jacobianMatrix[1, 1] += nurbs.NurbsDerivativeValuesHeta[k, j] * shellElement.ControlPoints[k].Y;
                jacobianMatrix[1, 2] += nurbs.NurbsDerivativeValuesHeta[k, j] * shellElement.ControlPoints[k].Z;
            }

            return jacobianMatrix;
        }

        private IList<GaussLegendrePoint3D> CreateElementGaussPoints(Element element)
        {
            GaussQuadrature gauss = new GaussQuadrature();
            return gauss.CalculateElementGaussPoints(element.Patch.DegreeKsi, element.Patch.DegreeHeta, element.Knots);
        }
        
        public IMatrix MassMatrix(IElement element)
        {
            throw new NotImplementedException();
        }

        public IMatrix DampingMatrix(IElement element)
        {
            throw new NotImplementedException();
        }

        public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
        {
            throw new NotImplementedException();
        }

        public IList<IList<IDofType>> GetElementDOFTypes(IElement element)
        {
            var nurbsElement = (NURBSKirchhoffLoveShellElement)element;
            dofTypes = new IDofType[nurbsElement.ControlPoints.Count][];
            for (int i = 0; i < nurbsElement.ControlPoints.Count; i++)
            {
                dofTypes[i] = controlPointDOFTypes;
            }
            return dofTypes;
        }

    }
}
