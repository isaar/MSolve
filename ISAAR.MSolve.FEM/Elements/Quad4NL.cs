using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Elements
{
    public class Quad4NL : IStructuralFiniteElement, IEmbeddedHostElement
    {
        protected static double determinantTolerance = 0.00000001;

        protected static int iInt = 2;
        protected static int iInt2 = iInt * iInt;
        protected static int iInt3 = iInt2 * iInt;
        protected readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y }; //Panos - 2 DoF per node
        protected readonly static DOFType[][] dofTypes = new DOFType[][] { nodalDOFTypes, nodalDOFTypes,
            nodalDOFTypes,nodalDOFTypes}; //Panos - 4 nodes only

        //protected readonly IFiniteElementMaterial3D[] materialsAtGaussPoints; //Panos - What about 2D Material? Implemented
        protected readonly ElasticMaterial2D[] materialsAtGaussPoints;

        protected IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();
        protected Quad4NL()
        {
        }
        public Quad4NL(IFiniteElementMaterial material)
        {
            materialsAtGaussPoints = new ElasticMaterial2D[iInt2];
            for (int i = 0; i < iInt2; i++)
                materialsAtGaussPoints[i] = (ElasticMaterial2D)material.Clone();
        }
        public Quad4NL(IFiniteElementMaterial material, IElementDOFEnumerator dofEnumerator)
            : this(material)
        {
            this.dofEnumerator = dofEnumerator;
        }
        public IElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }
        
        public double Density { get; set; }
        public double RayleighAlpha { get; set; }
        public double RayleighBeta { get; set; }
        
        protected double[,] GetCoordinates(IElement element)
        {
            //double[,] faXYZ = new double[dofTypes.Length, 3];
            double[,] faXY = new double[dofTypes.Length, 2];
            for (int i = 0; i < dofTypes.Length; i++)
            {
                faXY[i, 0] = element.INodes[i].X;
                faXY[i, 1] = element.INodes[i].Y;
            }
            return faXY;
        }
        protected double[,] GetCoordinatesTranspose(Element element)
        {
            double[,] faXY = new double[2, dofTypes.Length];
            for (int i = 0; i < dofTypes.Length; i++)
            {
                faXY[0, i] = element.Nodes[i].X;
                faXY[1, i] = element.Nodes[i].Y;
            }
            return faXY;
        }

        #region IElementType Members

        public int ID
        {
            //Change Element ID to random 14? I am not sure why is this needed
            get;
            set;
        }

        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.TwoD; }
        }

        public virtual IList<IList<DOFType>> GetElementDOFTypes(IElement element)
        {
            return dofTypes;
        }
        public IList<Node> GetNodesForMatrixAssembly(Element element)
        {
            return element.Nodes;
        }

        private double[] CalculateShapeFunctions(double coordinateKsi, double coordinateHeta)
        {
            // QUAD4 ELEMENT
            //
            //            Eta
            //            ^
            //#4(-1,1)    |     #3(1,1)
            //      ._____|_____.
            //      |     |     |
            //      |     |     |
            //      |     |_____|____>Xi (or psi - the same)
            //      |           |
            //      |           |
            //      .___________.
            // 
            //#1(-1,-1)        #2(1,-1)
            //
            double shapeFunctionKsi1 = (1.0 - coordinateKsi);
            double shapeFunctionKsi2 = (1.0 + coordinateKsi);
            double shapeFunctionHeta1 = (1.0 - coordinateHeta);
            double shapeFunctionHeta2 = (1.0 + coordinateHeta);
            
            var N1 = 0.25 * shapeFunctionKsi1 * shapeFunctionHeta1;
            var N2 = 0.25 * shapeFunctionKsi2 * shapeFunctionHeta1;
            var N3 = 0.25 * shapeFunctionKsi2 * shapeFunctionHeta2;
            var N4 = 0.25 * shapeFunctionKsi1 * shapeFunctionHeta2;
            return new[] {N1, N2, N3, N4};

        }

        private double[] CalculateShapeFunctionDerivatives(double coordinateKsi, double coordinateHeta)
        {
            double shapeFunctionKsi1 = (1.0 - coordinateKsi);
            double shapeFunctionKsi2 = (1.0 + coordinateKsi); 
            double shapeFunctionHeta1 = (1.0 - coordinateHeta);
            double shapeFunctionHeta2 = (1.0 + coordinateHeta);

            double[] shapeFunctionDerivatives = new double[8];

            shapeFunctionDerivatives[0] = -0.25*shapeFunctionHeta1; 
            shapeFunctionDerivatives[1] = 0.25* shapeFunctionHeta1; 
            shapeFunctionDerivatives[2] = 0.25* shapeFunctionHeta2;  
            shapeFunctionDerivatives[3] = -0.25* shapeFunctionHeta2; 

            shapeFunctionDerivatives[4] = -0.25* shapeFunctionKsi1;
            shapeFunctionDerivatives[5] = -0.25* shapeFunctionKsi2;
            shapeFunctionDerivatives[6] = 0.25* shapeFunctionKsi2;
            shapeFunctionDerivatives[7] = 0.25* shapeFunctionKsi1; 

            return shapeFunctionDerivatives;
        }

        private double[,] CalculateJacobianMatrix(double[,] nodeCoordinates, double[] shapeFunctionDerivatives)
        {
            double[,] jacobianMatrix = new double[2, 2];
            jacobianMatrix[0, 0] = shapeFunctionDerivatives[0] * nodeCoordinates[0, 0] + 
                                   shapeFunctionDerivatives[1] * nodeCoordinates[1, 0] + 
                                   shapeFunctionDerivatives[2] * nodeCoordinates[2, 0] + 
                                   shapeFunctionDerivatives[3] * nodeCoordinates[3, 0];
            jacobianMatrix[0, 1] = shapeFunctionDerivatives[0] * nodeCoordinates[0, 1] + 
                                   shapeFunctionDerivatives[1] * nodeCoordinates[1, 1] + 
                                   shapeFunctionDerivatives[2] * nodeCoordinates[2, 1] + 
                                   shapeFunctionDerivatives[3] * nodeCoordinates[3, 1];
            jacobianMatrix[1, 0] = shapeFunctionDerivatives[4] * nodeCoordinates[0, 0] + 
                                   shapeFunctionDerivatives[5] * nodeCoordinates[1, 0] +
                                   shapeFunctionDerivatives[6] * nodeCoordinates[2, 0] + 
                                   shapeFunctionDerivatives[7] * nodeCoordinates[3, 0];
            jacobianMatrix[1, 1] = shapeFunctionDerivatives[4] * nodeCoordinates[0, 1] + 
                                   shapeFunctionDerivatives[5] * nodeCoordinates[1, 1] + 
                                   shapeFunctionDerivatives[6] * nodeCoordinates[2, 1] + 
                                   shapeFunctionDerivatives[7] * nodeCoordinates[3, 1];

            return jacobianMatrix;
        }
        private double CalculateJacobianDeterminant(double[,] nodeCoordinates, double coordinateKsi, double coordinateHeta) 
        {
            double jacobianDeterminant =
                (nodeCoordinates[0, 0] * ((1 - coordinateHeta) * nodeCoordinates[1, 1] +
                                          (coordinateHeta - coordinateKsi) * nodeCoordinates[2, 1] +
                                          (coordinateKsi - 1) * nodeCoordinates[3, 1])) +
                (nodeCoordinates[1, 0] * ((coordinateHeta - 1) * nodeCoordinates[0, 1] +
                                          (coordinateKsi + 1) * nodeCoordinates[2, 1] +
                                          (-coordinateKsi - coordinateHeta) * nodeCoordinates[3, 1])) +
                (nodeCoordinates[2, 0] * ((coordinateKsi - coordinateHeta) * nodeCoordinates[0, 1] +
                                          (-coordinateKsi - 1) * nodeCoordinates[1, 1] +
                                          (coordinateHeta + 1) * nodeCoordinates[3, 1])) +
                (nodeCoordinates[3, 0] * ((1 - coordinateKsi) * nodeCoordinates[0, 1] +
                                          (coordinateKsi + coordinateHeta) * nodeCoordinates[1, 1] +
                                          (-coordinateHeta - 1) * nodeCoordinates[2, 1]));

            jacobianDeterminant = jacobianDeterminant*1/8;

            if (jacobianDeterminant < determinantTolerance)
                throw new ArgumentException(String.Format("Jacobian determinant is negative or under tolerance ({0} < {1})." +
                                                          "Check the order of nodes or the element geometry.", jacobianDeterminant, determinantTolerance));
            return jacobianDeterminant;
        }
        private double[,] CalculateJacobianInverseMatrix(double jacobianDeterminant, double[,] jacobianMatrix)
        {
            double fDetInv = 1.0 / jacobianDeterminant;

            double[,] jacobianInverseMatrix = new double[2, 2];
            jacobianInverseMatrix[0, 0] = (jacobianMatrix[1, 1]) / fDetInv;
            jacobianInverseMatrix[0, 1] = (-jacobianMatrix[0, 1]) / fDetInv;
            jacobianInverseMatrix[1, 0] = (-jacobianMatrix[1, 0]) / fDetInv;
            jacobianInverseMatrix[1, 1] = (jacobianMatrix[0, 0]) / fDetInv;

            return jacobianInverseMatrix;
        }

        private double[,] CalculateDeformationMatrix(double[,] nodeCoordinates, double[] shapeFunctionDerivatives, double coordinateKsi, double coordinateHeta, double jacobianDeterminant) //Panos - Calculate Deformation matrix B
        {
            double Aparameter;
            double Bparameter;
            double Cparameter;
            double Dparameter;

            Aparameter = 0.25 * (nodeCoordinates[0, 1] * (coordinateKsi - 1) + nodeCoordinates[1, 1] * (-1 - coordinateKsi) + nodeCoordinates[2, 1] * (1 + coordinateKsi) + nodeCoordinates[3, 1] * (1 - coordinateKsi));
            Bparameter = 0.25 * (nodeCoordinates[0, 1] * (coordinateHeta - 1) + nodeCoordinates[1, 1] * (1 - coordinateHeta) + nodeCoordinates[2, 1] * (1 + coordinateHeta) + nodeCoordinates[3, 1] * (-1 - coordinateHeta));
            Cparameter = 0.25 * (nodeCoordinates[0, 0] * (coordinateHeta - 1) + nodeCoordinates[1, 0] * (1 - coordinateHeta) + nodeCoordinates[2, 0] * (1 + coordinateHeta) + nodeCoordinates[3, 0] * (-1 - coordinateHeta));
            Dparameter = 0.25 * (nodeCoordinates[0, 0] * (coordinateKsi - 1) + nodeCoordinates[1, 0] * (-1 - coordinateKsi) + nodeCoordinates[2, 0] * (1 + coordinateKsi) + nodeCoordinates[3, 0] * (1 - coordinateKsi));

            double[,] Bmatrix = new double[3, 8];

            Bmatrix[0, 0] = (Aparameter * shapeFunctionDerivatives[0] - Bparameter * shapeFunctionDerivatives[4])/ jacobianDeterminant;
            Bmatrix[1, 0] = 0;
            Bmatrix[2, 0] = (Cparameter * shapeFunctionDerivatives[4] - Dparameter * shapeFunctionDerivatives[0])/ jacobianDeterminant;
            Bmatrix[0, 1] = 0;
            Bmatrix[1, 1] = (Cparameter * shapeFunctionDerivatives[4] - Dparameter * shapeFunctionDerivatives[0])/ jacobianDeterminant;
            Bmatrix[2, 1] = (Aparameter * shapeFunctionDerivatives[0] - Bparameter * shapeFunctionDerivatives[4])/ jacobianDeterminant;

            Bmatrix[0, 2] = (Aparameter * shapeFunctionDerivatives[1] - Bparameter * shapeFunctionDerivatives[5]) / jacobianDeterminant;
            Bmatrix[1, 2] = 0;
            Bmatrix[2, 2] = (Cparameter * shapeFunctionDerivatives[5] - Dparameter * shapeFunctionDerivatives[1]) / jacobianDeterminant;
            Bmatrix[0, 3] = 0;
            Bmatrix[1, 3] = (Cparameter * shapeFunctionDerivatives[5] - Dparameter * shapeFunctionDerivatives[1]) / jacobianDeterminant;
            Bmatrix[2, 3] = (Aparameter * shapeFunctionDerivatives[1] - Bparameter * shapeFunctionDerivatives[5]) / jacobianDeterminant;

            Bmatrix[0, 4] = (Aparameter * shapeFunctionDerivatives[2] - Bparameter * shapeFunctionDerivatives[6]) / jacobianDeterminant;
            Bmatrix[1, 4] = 0;
            Bmatrix[2, 4] = (Cparameter * shapeFunctionDerivatives[6] - Dparameter * shapeFunctionDerivatives[2]) / jacobianDeterminant;
            Bmatrix[0, 5] = 0;
            Bmatrix[1, 5] = (Cparameter * shapeFunctionDerivatives[6] - Dparameter * shapeFunctionDerivatives[2]) / jacobianDeterminant;
            Bmatrix[2, 5] = (Aparameter * shapeFunctionDerivatives[2] - Bparameter * shapeFunctionDerivatives[6]) / jacobianDeterminant;

            Bmatrix[0, 6] = (Aparameter * shapeFunctionDerivatives[3] - Bparameter * shapeFunctionDerivatives[7]) / jacobianDeterminant;
            Bmatrix[1, 6] = 0;
            Bmatrix[2, 6] = (Cparameter * shapeFunctionDerivatives[7] - Dparameter * shapeFunctionDerivatives[3]) / jacobianDeterminant;
            Bmatrix[0, 7] = 0;
            Bmatrix[1, 7] = (Cparameter * shapeFunctionDerivatives[7] - Dparameter * shapeFunctionDerivatives[3]) / jacobianDeterminant;
            Bmatrix[2, 7] = (Aparameter * shapeFunctionDerivatives[3] - Bparameter * shapeFunctionDerivatives[7]) / jacobianDeterminant;

            return Bmatrix;
        }

        private GaussLegendrePoint3D[] CalculateGaussMatrices(double[,] nodeCoordinates)
        {
            GaussLegendrePoint1D[] integrationPointsPerAxis =
                GaussQuadrature.GetGaussLegendrePoints(iInt);
            int totalSamplingPoints = (int)Math.Pow(integrationPointsPerAxis.Length, 2);

            GaussLegendrePoint3D[] integrationPoints = new GaussLegendrePoint3D[totalSamplingPoints];

            int counter = -1;
            foreach (GaussLegendrePoint1D pointXi in integrationPointsPerAxis)
            {
                foreach (GaussLegendrePoint1D pointEta in integrationPointsPerAxis)
                {
                    counter += 1;
                    double ksi = pointXi.Coordinate;
                    double heta = pointEta.Coordinate;
                    double[] shapeFunctionDerivatives = this.CalculateShapeFunctionDerivatives(ksi, heta);
                    double fDetJ = this.CalculateJacobianDeterminant(nodeCoordinates, ksi, heta);
                    double[,] deformationMatrix = this.CalculateDeformationMatrix(nodeCoordinates, shapeFunctionDerivatives, ksi, heta, fDetJ);
                    double weightFactor = pointXi.WeightFactor * pointEta.WeightFactor * fDetJ; 
                    integrationPoints[counter] = new GaussLegendrePoint3D(ksi, heta,0, deformationMatrix, weightFactor);
                }
            }
            return integrationPoints;
        }

        public virtual IMatrix2D StiffnessMatrix(IElement element)
        {
            double[,] coordinates = this.GetCoordinates(element);
            GaussLegendrePoint3D[] integrationPoints = this.CalculateGaussMatrices(coordinates);
            SymmetricMatrix2D stiffnessMatrix = new SymmetricMatrix2D(8);
            int pointId = -1;
            foreach (GaussLegendrePoint3D gaussPoint in integrationPoints)
            {
                pointId++;
                IMatrix2D constitutiveMatrix = materialsAtGaussPoints[pointId].ConstitutiveMatrix;
                double[,] b = gaussPoint.DeformationMatrix;
                for (int i = 0; i < 8; i++)
                {
                    double[] eb = new double[3];
                    for (int iE = 0; iE < 3; iE++)
                    {
                        eb[iE] = (constitutiveMatrix[iE, 0] * b[0, i]) + (constitutiveMatrix[iE, 1] * b[1, i]) +
                                 (constitutiveMatrix[iE, 2] * b[2, i]);
                    }


                    for (int j = i; j < 8; j++)
                    {
                        double stiffness = (b[0, j] * eb[0]) + (b[1, j] * eb[1]) + (b[2, j] * eb[2]);
                        stiffnessMatrix[i, j] += stiffness * gaussPoint.WeightFactor;
                    }
                }
            }

            return stiffnessMatrix;
        }
        
        #endregion

        public virtual IMatrix2D MassMatrix(IElement element)
        {
            return CalculateConsistentMass(element);
        }

        public virtual IMatrix2D DampingMatrix(IElement element)
        {
            var m = MassMatrix(element);
            var lc = m as ILinearlyCombinable;
            lc.LinearCombination(new double[] { RayleighAlpha, RayleighBeta }, new IMatrix2D[] { MassMatrix(element), StiffnessMatrix(element) });
            return m;
        }

        public IMatrix2D CalculateConsistentMass(IElement element)
        {
            double[,] coordinates = this.GetCoordinates(element);
            GaussLegendrePoint3D[] integrationPoints = this.CalculateGaussMatrices(coordinates);
            SymmetricMatrix2D consistentMass = new SymmetricMatrix2D(8);
            foreach (GaussLegendrePoint3D gaussPoint in integrationPoints)
            {
                double[] shapeFunctionValues = CalculateShapeFunctions(gaussPoint.Xi, gaussPoint.Eta);
                double weightDensity = gaussPoint.WeightFactor * Density;
                for (int iShapeFunction = 0; iShapeFunction < shapeFunctionValues.Length; iShapeFunction++)
                {
                    for (int jShapeFunction = iShapeFunction; jShapeFunction < shapeFunctionValues.Length; jShapeFunction++)
                    {
                        consistentMass[2 * iShapeFunction, 2 * jShapeFunction]
                            += shapeFunctionValues[iShapeFunction] *
                              shapeFunctionValues[jShapeFunction] *
                              weightDensity;
                    }
                    for (int jShapeFunction = iShapeFunction; jShapeFunction < shapeFunctionValues.Length; jShapeFunction++)
                    {
                        consistentMass[(2 * iShapeFunction) + 1, (2 * jShapeFunction) + 1] =
                            consistentMass[2 * iShapeFunction, 2 * jShapeFunction];
                    }
                }
            }
            return consistentMass;
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            double[,] faXY = GetCoordinates(element);
            double[,] faDS = new double[iInt3, 24];
            double[,] faS = new double[iInt3, 8];
            double[,,] faB = new double[iInt3, 24, 6];
            double[] faDetJ = new double[iInt3];
            double[,,] faJ = new double[iInt3, 3, 3];
            double[] faWeight = new double[iInt3];
            double[,] fadStrains = new double[iInt3, 6];
            double[,] faStrains = new double[iInt3, 6];
            //	CalcQ4GaussMatrices(ref iInt, nodeCoordinates, faWeight, faS, shapeFunctionDerivatives, jacobianMatrix, faDetJ, faB);
            //	CalcQ4Strains(ref iInt, faB, localDisplacements, faStrains);
            //	CalcQ4Strains(ref iInt, faB, localdDisplacements, fadStrains);

            double[] dStrains = new double[6];
            double[] strains = new double[6];
            for (int i = 0; i < materialsAtGaussPoints.Length; i++)
            {
                for (int j = 0; j < 6; j++) dStrains[j] = fadStrains[i, j];
                for (int j = 0; j < 6; j++) strains[j] = faStrains[i, j];
                materialsAtGaussPoints[i].UpdateMaterial(new StressStrainVectorContinuum2D(dStrains));
            }

            return new Tuple<double[], double[]>(strains, materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses.Data);
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public double[] CalculateForces(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            //Vector<double> d = new Vector<double>(localdDisplacements.Length);
            //for (int i = 0; i < localdDisplacements.Length; i++) 
            //    //d[i] = localdDisplacements[i] + localTotalDisplacements[i];
            //    d[i] = localTotalDisplacements[i];
            //double[] faForces = new double[24];
            //StiffnessMatrix(element).Multiply(d, faForces);

            double[,] faStresses = new double[iInt3, 6];
            for (int i = 0; i < materialsAtGaussPoints.Length; i++)
                for (int j = 0; j < 6; j++) faStresses[i, j] = materialsAtGaussPoints[i].Stresses[j];

            double[,] faXYZ = GetCoordinates(element);
            double[,] faDS = new double[iInt3, 24];
            double[,] faS = new double[iInt3, 8];
            double[,,] faB = new double[iInt3, 24, 6];
            double[] faDetJ = new double[iInt3];
            double[,,] faJ = new double[iInt3, 3, 3];
            double[] faWeight = new double[iInt3];
            double[] faForces = new double[24];
            //	CalcQ4GaussMatrices(ref iInt, faXYZ, faWeight, faS, shapeFunctionDerivatives, jacobianMatrix, faDetJ, faB);
            //	CalcQ4Forces(ref iInt, faB, faWeight, faStresses, faForces);

            return faForces;
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            Vector accelerations = new Vector(8);
            IMatrix2D massMatrix = MassMatrix(element);

            foreach (MassAccelerationLoad load in loads)
            {
                int index = 0;
                foreach (DOFType[] nodalDOFTypes in dofTypes)
                    foreach (DOFType dofType in nodalDOFTypes)
                    {
                        if (dofType == load.DOF) accelerations[index] += load.Amount;
                        index++;
                    }
            }
            double[] forces = new double[8];
            massMatrix.Multiply(accelerations, forces);
            return forces;
        }

        public void ClearMaterialState()
        {
            throw new NotImplementedException();
        }

        public void SaveMaterialState()
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialStresses()
        {
            throw new NotImplementedException();
        }

        #region IStructuralFiniteElement Members

        public bool MaterialModified
        {
            get
            {
                throw new NotImplementedException();
            }
        }

        public void ResetMaterialModified()
        {
            throw new NotImplementedException();
        }

        #endregion

        #region IEmbeddedHostElement Members

        public EmbeddedNode BuildHostElementEmbeddedNode(Element element, Node node, IEmbeddedDOFInHostTransformationVector transformationVector)
        {
            var points = GetNaturalCoordinates(element, node);
            if (points.Length == 0) return null;

            element.EmbeddedNodes.Add(node);
            var embeddedNode = new EmbeddedNode(node, element, transformationVector.GetDependentDOFTypes);
            for (int i = 0; i < points.Length; i++)
                embeddedNode.Coordinates.Add(points[i]);
            return embeddedNode;
        }

        public double[] GetShapeFunctionsForNode(Element element, EmbeddedNode node)
        {
            double[,] elementCoordinates = GetCoordinatesTranspose(element);
            var shapeFunctions = CalculateShapeFunctions(node.Coordinates[0], node.Coordinates[1]);
            var ShapeFunctionDerivatives = CalculateShapeFunctionDerivatives(node.Coordinates[0], node.Coordinates[1]);
            var jacobian = CalculateJacobianMatrix(elementCoordinates, ShapeFunctionDerivatives);

            return new double[]
            {
                shapeFunctions[0], shapeFunctions[1], shapeFunctions[2], shapeFunctions[3],
                ShapeFunctionDerivatives[0], ShapeFunctionDerivatives[1], ShapeFunctionDerivatives[2], ShapeFunctionDerivatives[3], ShapeFunctionDerivatives[4], ShapeFunctionDerivatives[5], ShapeFunctionDerivatives[6], ShapeFunctionDerivatives[7],
                jacobian[0, 0], jacobian[0, 1], jacobian[1, 0], jacobian[1, 1]
            };

            //double[] coords = GetNaturalCoordinates(element, node);
            //return CalcH8Shape(coords[0], coords[1], coords[2]);

            //double fXiP = (1.0 + coords[0]) * 0.5;
            //double fEtaP = (1.0 + coords[1]) * 0.5;
            //double fZetaP = (1.0 + coords[2]) * 0.5;
            //double fXiM = (1.0 - coords[0]) * 0.5;
            //double fEtaM = (1.0 - coords[1]) * 0.5;
            //double fZetaM = (1.0 - coords[2]) * 0.5;

            //return new double[] { fXiM * fEtaM * fZetaM,
            //    fXiP * fEtaM * fZetaM,
            //    fXiP * fEtaP * fZetaM,
            //    fXiM * fEtaP * fZetaM,
            //    fXiM * fEtaM * fZetaP,
            //    fXiP * fEtaM * fZetaP,
            //    fXiP * fEtaP * fZetaP,
            //    fXiM * fEtaP * fZetaP };
        }

        private double[] GetNaturalCoordinates(Element element, Node node)
        {
            double[] mins = new double[] { element.Nodes[0].X, element.Nodes[0].Y };
            double[] maxes = new double[] { element.Nodes[0].X, element.Nodes[0].Y };
            for (int i = 0; i < element.Nodes.Count; i++)
            {
                mins[0] = mins[0] > element.Nodes[i].X ? element.Nodes[i].X : mins[0];
                mins[1] = mins[1] > element.Nodes[i].Y ? element.Nodes[i].Y : mins[1];
                maxes[0] = maxes[0] < element.Nodes[i].X ? element.Nodes[i].X : maxes[0];
                maxes[1] = maxes[1] < element.Nodes[i].Y ? element.Nodes[i].Y : maxes[1];
            }
            //return new double[] { (node.X - mins[0]) / ((maxes[0] - mins[0]) / 2) - 1,
            //    (node.Y - mins[1]) / ((maxes[1] - mins[1]) / 2) - 1,
            //    (node.Z - mins[2]) / ((maxes[2] - mins[2]) / 2) - 1 };

            bool maybeInsideElement = node.X <= maxes[0] && node.X >= mins[0] &&
                                          node.Y <= maxes[1] && node.Y >= mins[1];
            if (maybeInsideElement == false) return new double[0];


            const int jacobianSize = 3;
            const int maxIterations = 1000;
            const double tolerance = 1e-10;
            int iterations = 0;
            double deltaNaturalCoordinatesNormSquare = 100;
            double[] naturalCoordinates = new double[] { 0, 0, 0 };
            const double toleranceSquare = tolerance * tolerance;

            while (deltaNaturalCoordinatesNormSquare > toleranceSquare && iterations < maxIterations)
            {
                iterations++;
                var shapeFunctions = CalculateShapeFunctions(naturalCoordinates[0], naturalCoordinates[1]);
                double[] coordinateDifferences = new double[] { 0, 0 };
                for (int i = 0; i < shapeFunctions.Length; i++)
                {
                    coordinateDifferences[0] += shapeFunctions[i] * element.Nodes[i].X;
                    coordinateDifferences[1] += shapeFunctions[i] * element.Nodes[i].Y;
                }
                coordinateDifferences[0] = node.X - coordinateDifferences[0];
                coordinateDifferences[1] = node.Y - coordinateDifferences[1];

                double[,] faXY = GetCoordinatesTranspose(element);
                double[] ShapeFunctionDerivatives = CalculateShapeFunctionDerivatives(naturalCoordinates[0], naturalCoordinates[1]);
                //SOS PANOS - THIS IS WRONG, I SHOULD IMPLEMENT INVERSE JACOBIAN
                var inverseJacobian = CalculateJacobianMatrix(faXY, ShapeFunctionDerivatives);
                //SOS PANOS END
                double[] deltaNaturalCoordinates = new double[] { 0, 0, 0 };
                for (int i = 0; i < jacobianSize; i++)
                    for (int j = 0; j < jacobianSize; j++)
                        deltaNaturalCoordinates[i] += inverseJacobian[j, i] * coordinateDifferences[j];
                for (int i = 0; i < 3; i++)
                    naturalCoordinates[i] += deltaNaturalCoordinates[i];

                deltaNaturalCoordinatesNormSquare = 0;
                for (int i = 0; i < 3; i++)
                    deltaNaturalCoordinatesNormSquare += deltaNaturalCoordinates[i] * deltaNaturalCoordinates[i];
                //deltaNaturalCoordinatesNormSquare = Math.Sqrt(deltaNaturalCoordinatesNormSquare);
            }

            return naturalCoordinates.Count(x => Math.Abs(x) - 1.0 > tolerance) > 0 ? new double[0] : naturalCoordinates;
        }

        #endregion
    }
}
