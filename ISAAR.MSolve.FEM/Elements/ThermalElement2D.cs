using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.FEM.Interpolation.Jacobians;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

//TODO: Is there any point in having different material properties per Gauss point?
namespace ISAAR.MSolve.FEM.Elements
{
    public class ThermalElement2D : IFiniteElement, IEmbeddedHostElement
    {
        private readonly DOFType[][] dofTypes; //TODO: this should not be stored for each element. Instead store it once for each Quad4, Tri3, etc. Otherwise create it on the fly.
        private readonly ThermalMaterial material;
        //private readonly Dictionary<GaussPoint2D, ThermalMaterial> materialsAtGaussPoints;


        public ThermalElement2D(double thickness, IReadOnlyList<Node2D> nodes, IIsoparametricInterpolation2D interpolation,
            IQuadrature2D quadratureForStiffness, IQuadrature2D quadratureForConsistentMass,
            IGaussPointExtrapolation2D gaussPointExtrapolation,
            ThermalMaterial material)
        {
            this.material = material;
            this.GaussPointExtrapolation = gaussPointExtrapolation;
            this.Nodes = nodes;
            this.Interpolation = interpolation;
            this.QuadratureForConsistentMass = quadratureForConsistentMass;
            this.QuadratureForStiffness = quadratureForStiffness;
            this.Thickness = thickness;

            dofTypes = new DOFType[nodes.Count][];
            for (int i = 0; i < interpolation.NumFunctions; ++i) dofTypes[i] = new DOFType[] { DOFType.Temperature };
        }

        public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

        public int ID => throw new NotImplementedException(
            "Element type codes should be in a settings class. Even then it's a bad design choice");

        public IGaussPointExtrapolation2D GaussPointExtrapolation { get; }
        public IIsoparametricInterpolation2D Interpolation { get; }
        public IReadOnlyList<Node2D> Nodes { get; }
        public IQuadrature2D QuadratureForConsistentMass { get; }
        public IQuadrature2D QuadratureForStiffness { get; }
        public double Thickness { get; }

        public bool MaterialModified => throw new NotImplementedException();

        public IElementDOFEnumerator DOFEnumerator { get; set; } = new GenericDOFEnumerator();

        public IMatrix2D MassMatrix(IElement element)
        {
            return BuildCapacityMatrix();
        }

        public Matrix2D BuildCapacityMatrix()
        {
            int numDofs = Nodes.Count;
            var capacity = new Matrix2D(numDofs, numDofs);
            IReadOnlyList<Vector> shapeFunctions =
                Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
            IReadOnlyList<Matrix2D> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

            for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
            {
                Matrix2D shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
                Matrix2D partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
                var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
                double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
                capacity.AxpyIntoThis(partial, dA);
            }

            //WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            capacity.Scale(Thickness * material.Density * material.SpecialHeatCoeff);
            return capacity;
        }

        public Matrix2D BuildConductivityMatrix()
        {
            int numDofs = Nodes.Count;
            var conductivity = new Matrix2D(numDofs, numDofs);
            IReadOnlyList<Matrix2D> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

            for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
            {
                // Calculate the necessary quantities for the integration
                //Matrix2D constitutive = (Matrix2D)(materialsAtGaussPoints[gp].ConstitutiveMatrix); // ugly cast will be removed along with the legacy Matrix classes
                var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
                Matrix2D shapeGradientsCartesian =
                    jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
                Matrix2D deformation = BuildDeformationMatrix(shapeGradientsCartesian);

                // Contribution of this gauss point to the element stiffness matrix
                Matrix2D partialK = deformation.Transpose() * deformation;
                //Matrix2D partialΚ = deformation.Transpose() * (constitutive * deformation);
                //partialK.Scale(materialsAtGaussPoints[gaussPoint].ThermalConductivity);

                double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight; //TODO: this is used by all methods that integrate. I should cache it.
                conductivity.AxpyIntoThis(partialK, dA * material.ThermalConductivity);
            }

            conductivity.Scale(Thickness);
            return conductivity;
        }

        // Provatidis uses two distinct vectors K = N,x^T * k * N,x + N,y^T * k * N,y
        //private (Matrix2D dNdX, Matrix2D dNdY) CalcdNdx(EvalShapeGradients2D shapeGrad)
        //{
        //    int n = Nodes.Count;
        //    var dNdX = new double[n, 1];
        //    var dNdY = new double[n, 1];
        //    for (int i = 0; i < n; ++i)
        //    {
        //        dNdX[i, 0] = shapeGrad[i][0];
        //        dNdY[i, 0] = shapeGrad[i][1];
        //    }
        //    return (new Matrix2D(dNdX), new Matrix2D(dNdY));
        //}

        private Matrix2D BuildDeformationMatrix(Matrix2D shapeGradientsCartesian) 
        {
            //TODO: isn't this just the transpose of [dNi/dxj]?
            var deformation = new Matrix2D(2, Nodes.Count);
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                deformation[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 0];
                deformation[1, nodeIdx] = shapeGradientsCartesian[nodeIdx, 1];
            }
            return deformation;
        }

        /// <summary>
        /// The shape function matrix is 1-by-n, where n = is the number of shape functions.
        /// </summary>
        public Matrix2D BuildShapeFunctionMatrix(Vector shapeFunctions) //TODO: reconsider this. As it is, it just returns the shape functions in a Matrix2D
        {
            var array2D = new double[1, shapeFunctions.Length];
            for (int i = 0; i < shapeFunctions.Length; ++i)
            {
                array2D[0, i] = shapeFunctions[i];
            }
            return new Matrix2D(array2D);
        }

        public IList<IList<DOFType>> GetElementDOFTypes(IElement element) => dofTypes;

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

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            throw new NotImplementedException();
        }

        public void SaveMaterialState()
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialState()
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialStresses()
        {
            throw new NotImplementedException();
        }

        public IMatrix2D StiffnessMatrix(IElement element)
        {
            return DOFEnumerator.GetTransformedMatrix(BuildConductivityMatrix());
        }

        public IMatrix2D DampingMatrix(IElement element)
        {
            throw new NotImplementedException();
        }

        public EmbeddedNode BuildHostElementEmbeddedNode(Element element, Node node, IEmbeddedDOFInHostTransformationVector transformationVector)
        {
            IInverseInterpolation2D inverseInterpolation = Interpolation.CreateInverseMappingFor(Nodes);
            double[] naturalCoordinates = inverseInterpolation.TransformPointCartesianToNatural(new CartesianPoint2D(node.X, node.Y)).Coordinates;

            if (naturalCoordinates.Length == 0) return null;

            element.EmbeddedNodes.Add(node);
            var embeddedNode = new EmbeddedNode(node, element, transformationVector.GetDependentDOFTypes);
            for (int i = 0; i < naturalCoordinates.Length; i++)
                embeddedNode.Coordinates.Add(naturalCoordinates[i]);
            return embeddedNode;
        }

        public double[] GetShapeFunctionsForNode(Element element, EmbeddedNode node)
        {
            return Interpolation.EvaluateFunctionsAt(new NaturalPoint2D(node.Coordinates[0], node.Coordinates[1])).Data;

            //TODO: This method originally returned an array containing the shape functions, shape function derivatives and entries of the inverse jacobian matrix.
            //      a) This is retarded and extremely difficult to work with. b) This array is used in Hexa8TranslationAndRotationTransformationVector. Cohesive embedding 
            //      and regular embedding in thermal elements only need the shape functions, so we will return these for now.

            //double[,] elementCoordinates = GetCoordinatesTranspose(element);
            //var shapeFunctions = CalcH8Shape(node.Coordinates[0], node.Coordinates[1], node.Coordinates[2]);
            //var nablaShapeFunctions = CalcH8NablaShape(node.Coordinates[0], node.Coordinates[1], node.Coordinates[2]);
            //var jacobian = CalcH8JDetJ(elementCoordinates, nablaShapeFunctions);

            //return new double[]
            //{
            //    shapeFunctions[0], shapeFunctions[1], shapeFunctions[2], shapeFunctions[3], shapeFunctions[4], shapeFunctions[5], shapeFunctions[6], shapeFunctions[7],
            //    nablaShapeFunctions[0], nablaShapeFunctions[1], nablaShapeFunctions[2], nablaShapeFunctions[3], nablaShapeFunctions[4], nablaShapeFunctions[5], nablaShapeFunctions[6], nablaShapeFunctions[7],
            //    nablaShapeFunctions[8], nablaShapeFunctions[9], nablaShapeFunctions[10], nablaShapeFunctions[11], nablaShapeFunctions[12], nablaShapeFunctions[13], nablaShapeFunctions[14], nablaShapeFunctions[15],
            //    nablaShapeFunctions[16], nablaShapeFunctions[17], nablaShapeFunctions[18], nablaShapeFunctions[19], nablaShapeFunctions[20], nablaShapeFunctions[21], nablaShapeFunctions[22], nablaShapeFunctions[23],
            //    jacobian.Item1[0, 0], jacobian.Item1[0, 1], jacobian.Item1[0, 2], jacobian.Item1[1, 0], jacobian.Item1[1, 1], jacobian.Item1[1, 2], jacobian.Item1[2, 0], jacobian.Item1[2, 1], jacobian.Item1[2, 2],
            //    jacobian.Item2[0, 0], jacobian.Item2[0, 1], jacobian.Item2[0, 2], jacobian.Item2[1, 0], jacobian.Item2[1, 1], jacobian.Item2[1, 2], jacobian.Item2[2, 0], jacobian.Item2[2, 1], jacobian.Item2[2, 2]
            //};
        }
    }
}
