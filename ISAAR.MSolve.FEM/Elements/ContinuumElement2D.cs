using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

//TODO: Damping matrix calculation needs redesign for all of MSolve. For this class, see DampingMatrix().
//TODO: Materials also need redesign. Some properties are the same for all instances of a material class, some are the same for
//      all Gauss points of an element but differ across elements, some differ per Gauss point but are constant, and others
//      differ per Gauss point and per analysis iteration.
//TODO: Why does http://kis.tu.kielce.pl/mo/COLORADO_FEM/colorado/IFEM.Ch31.pdf Fig 31.9 take half the thickness when computing
//      the consistent mass matrix of Quad4 elements? For Tri3, the full thickness is used, as seen in Fig 31.7
//TODO: Simple Tri3 elements are more efficient than isoparamateric Tri3 elements. Many optimizations could be made. See
//      https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch24.d/IFEM.Ch24.pdf for more.
//TODO: The shape of finite elements needs to be checked before the analysis: wrong node order, clockwise node order, the  
//      element's shape is too distorted, midpoints are too close to corners in quadratic elements, etc.
namespace ISAAR.MSolve.FEM.Elements
{
    /// <summary>
    /// Represents a continuum finite element for 2D problems. Specific elements (e.g. Quad4, Tri6, ...) can be created using
    /// the appropriate <see cref="IIsoparametricInterpolation2D"/>, <see cref="IQuadrature2D"/> etc. strategies. The thickness
    /// of this element is uniform, therefore it is necessary to use finer meshes to simulate domains with variable thickness.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ContinuumElement2D: IStructuralFiniteElement
    {
        private readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y };
        private readonly DOFType[][] dofTypes; //TODO: this should not be stored for each element. Instead store it once for each Quad4, Tri3, etc. Otherwise create it on the fly.
        private DynamicMaterial dynamicProperties;
        private readonly Dictionary<GaussPoint2D, ElasticMaterial2D> materialsAtGaussPoints;

        public ContinuumElement2D(double thickness, IReadOnlyList<Node2D> nodes, IIsoparametricInterpolation2D interpolation,
            IQuadrature2D quadratureForStiffness, IQuadrature2D quadratureForConsistentMass, 
            IGaussPointExtrapolation2D gaussPointExtrapolation, 
            Dictionary<GaussPoint2D, ElasticMaterial2D> materialsAtGaussPoints, DynamicMaterial dynamicProperties)
        {
            this.dynamicProperties = dynamicProperties;
            this.materialsAtGaussPoints = materialsAtGaussPoints;
            this.GaussPointExtrapolation = gaussPointExtrapolation;
            this.Nodes = nodes;
            this.Interpolation = interpolation;
            this.QuadratureForConsistentMass = quadratureForConsistentMass;
            this.QuadratureForStiffness = quadratureForStiffness;
            this.Thickness = thickness;

            dofTypes = new DOFType[nodes.Count][];
            for (int i = 0; i < interpolation.NumFunctions; ++i) dofTypes[i] = new DOFType[] { DOFType.X, DOFType.Y };
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

        public Matrix2D BuildConsistentMassMatrix()
        {
            int numDofs = 2 * Nodes.Count;
            var mass = new Matrix2D(numDofs, numDofs);
            Dictionary<GaussPoint2D, EvalInterpolation2D> evalInterpolations =
                Interpolation.EvaluateAllAtGaussPoints(Nodes, QuadratureForConsistentMass);

            foreach (GaussPoint2D gaussPoint in QuadratureForConsistentMass.IntegrationPoints)
            {
                Matrix2D shapeFunctionMatrix = evalInterpolations[gaussPoint].BuildShapeFunctionMatrix();
                Matrix2D partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
                double dA = evalInterpolations[gaussPoint].Jacobian.Determinant * gaussPoint.Weight;
                mass.AxpyIntoThis(partial, dA);
            }

            //WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            mass.Scale(Thickness * dynamicProperties.Density); 
            return mass;
        }

        public Matrix2D BuildLumpedMassMatrix()
        {
            int numDofs = 2 * Nodes.Count;
            var lumpedMass = new Matrix2D(numDofs, numDofs);
            Dictionary<GaussPoint2D, EvalShapeGradients2D> shapeGradients =
                Interpolation.EvaluateGradientsAtGaussPoints(Nodes, QuadratureForConsistentMass);

            // Contribution of each Gauss point to the element's area
            //TODO: Perhaps I could calculate the volume of the element without going through each Gauss point. Probably the 
            //      nodes are needed instead of the GPs. For linear elements I can find the area geometrically (as polygons).
            double area = 0;
            foreach (GaussPoint2D gaussPoint in QuadratureForConsistentMass.IntegrationPoints)
            {
                area += shapeGradients[gaussPoint].Jacobian.Determinant * gaussPoint.Weight;
            }

            // Divide the total mass uniformly for each node
            double nodalMass = Thickness * area * dynamicProperties.Density / Nodes.Count;
            for (int i = 0; i < numDofs; ++i) lumpedMass[i, i] = nodalMass;

            return lumpedMass;
        }

        public Matrix2D BuildStiffnessMatrix()
        {
            int numDofs = 2 * Nodes.Count;
            var stiffness = new Matrix2D(numDofs, numDofs);
            Dictionary<GaussPoint2D, EvalShapeGradients2D> shapeGradients =
                Interpolation.EvaluateGradientsAtGaussPoints(Nodes, QuadratureForStiffness);

            foreach (GaussPoint2D gaussPoint in QuadratureForStiffness.IntegrationPoints)
            {
                // Calculate the necessary quantities for the integration
                Matrix2D constitutive = (Matrix2D)(materialsAtGaussPoints[gaussPoint].ConstitutiveMatrix); // ugly cast will be removed along with the retarded legacy Matrix classes
                Matrix2D deformation = BuildDeformationMatrix(shapeGradients[gaussPoint]);

                // Contribution of this gauss point to the element stiffness matrix
                Matrix2D partial = deformation.Transpose() * (constitutive * deformation);
                double dA = shapeGradients[gaussPoint].Jacobian.Determinant * gaussPoint.Weight;
                stiffness.AxpyIntoThis(partial, dA);
            }
            stiffness.Scale(Thickness);
            return stiffness;
        }

        //TODO: I think this method must be removed from IFiniteElement altogether. This procedure shoud be done for the global 
        //      mass matrix, once at the start of the dynamic analysis. The resulting vectors for each direction of the ground 
        //      motion should be stored. Then at each timestep they only need to be scaled and added to the rhs vector. The mass 
        //      matrix doesn't change, so there is not reason to recompute it at each time step.
        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            int numDofs = 2 * Nodes.Count;
            Vector accelerations = new Vector(numDofs);
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
            double[] forces = new double[numDofs];
            massMatrix.Multiply(accelerations, forces);
            return forces;
        }

        public double[] CalculateForces(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        //TODO: this method is probably not necessary at all
        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        //TODO: this method must be changed. It should calculates strains, stresses at GPs or nodes.
        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, 
            double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialState()
        {
            foreach (ElasticMaterial2D m in materialsAtGaussPoints.Values) m.ClearState();
        }

        public void ClearMaterialStresses()
        {
            foreach (ElasticMaterial2D m in materialsAtGaussPoints.Values) m.ClearStresses();
        }

        public IMatrix2D DampingMatrix(IElement element)
        {
            //TODO: Stiffness and mass matrices have already been computed probably. Reuse them.
            //TODO: Perhaps with Rayleigh damping, the global damping matrix should be created directly from global mass and stiffness matrices.
            Matrix2D damping = BuildStiffnessMatrix();
            damping.Scale(dynamicProperties.RayleighCoeffStiffness);
            damping.AxpyIntoThis(MassMatrix(element), dynamicProperties.RayleighCoeffMass);
            return damping;
        }

        public IElementDOFEnumerator DOFEnumerator { get; set; } = new GenericDOFEnumerator();

        public IList<IList<DOFType>> GetElementDOFTypes(IElement element) => dofTypes;

        public IMatrix2D MassMatrix(IElement element)
        {
            return BuildConsistentMassMatrix();
            //return BuildLumpedMassMatrix();
        }

        public bool MaterialModified
        {
            get
            {
                foreach (ElasticMaterial2D material in materialsAtGaussPoints.Values)
                    if (material.Modified) return true;
                return false;
            }
        }

        public void ResetMaterialModified()
        {
            foreach (ElasticMaterial2D material in materialsAtGaussPoints.Values) material.ResetModified();
        }

        public void SaveMaterialState()
        {
            foreach (ElasticMaterial2D m in materialsAtGaussPoints.Values) m.SaveState();
        }

        //TODO: why do I need the wrapping element?
        public IMatrix2D StiffnessMatrix(IElement element)
        {
            return BuildStiffnessMatrix();
        }

        /// <summary>
        /// Calculate strains (exx, eyy, 2exy) and stresses (sxx, syy, sxy) at integration points, store them in the materials 
        /// and return them (e.g. for postprocessing). The order of the tensors is the same as the order of the integration 
        /// points defined by <see cref="QuadratureForStiffness"/>.
        /// </summary>
        /// <param name="localDisplacements"></param>
        /// <returns></returns>
        public (IReadOnlyList<double[]> strains, IReadOnlyList<double[]> stresses) UpdateStrainsStressesAtGaussPoints(
            double[] localDisplacements)
        {
            var localDisplVector = new Vector(localDisplacements);
            int numGPs = QuadratureForStiffness.IntegrationPoints.Count;
            var strains = new double[numGPs][];
            var stresses = new double[numGPs][];
            Dictionary<GaussPoint2D, EvalShapeGradients2D> shapeGradients =
                Interpolation.EvaluateGradientsAtGaussPoints(Nodes, QuadratureForStiffness);
            for (int i = 0; i < numGPs; ++i)
            {
                GaussPoint2D gaussPoint = QuadratureForStiffness.IntegrationPoints[i];
                IMatrix2D constitutive = materialsAtGaussPoints[gaussPoint].ConstitutiveMatrix;
                Matrix2D deformation = BuildDeformationMatrix(shapeGradients[gaussPoint]);
                strains[i] = new double[3];
                deformation.Multiply(localDisplVector, strains[i]);
                stresses[i] = new double[3];
                constitutive.Multiply(new Vector(strains[i]), stresses[i]);
            }
            return (strains, stresses);
        }

        private Matrix2D BuildDeformationMatrix(EvalShapeGradients2D shapeGradients)
        {
            var deformation = new Matrix2D(3, 2 * Nodes.Count);
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                int col1 = 2 * nodeIdx;
                int col2 = 2 * nodeIdx + 1;
                IReadOnlyList<double> dNdX = shapeGradients[nodeIdx];

                deformation[0, col1] = dNdX[0];
                deformation[1, col2] = dNdX[1];
                deformation[2, col1] = dNdX[1];
                deformation[2, col2] = dNdX[0];
            }
            return deformation;
        }
    }
}
