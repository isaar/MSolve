using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.FEM.Interpolation.Jacobians;
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
        private readonly IReadOnlyList<ElasticMaterial2D> materialsAtGaussPoints;

        public ContinuumElement2D(double thickness, IReadOnlyList<Node2D> nodes, IIsoparametricInterpolation2D interpolation,
            IQuadrature2D quadratureForStiffness, IQuadrature2D quadratureForConsistentMass, 
            IGaussPointExtrapolation2D gaussPointExtrapolation,
            IReadOnlyList<ElasticMaterial2D> materialsAtGaussPoints, DynamicMaterial dynamicProperties)
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
            IReadOnlyList<Matrix2D> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

            // Contribution of each Gauss point to the element's area
            //TODO: Perhaps I could calculate the volume of the element without going through each Gauss point. Probably the 
            //      nodes are needed instead of the GPs. For linear elements I can find the area geometrically (as polygons).
            //TODO: this should have been cached when integrating other quantities (e.g. stiffness)
            double area = 0;
            for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
            {
                var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
                area += jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
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
            IReadOnlyList<Matrix2D> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

            for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
            {
                // Calculate the necessary quantities for the integration
                Matrix2D constitutive = (Matrix2D)(materialsAtGaussPoints[gp].ConstitutiveMatrix); // ugly cast will be removed along with the retarded legacy Matrix classes
                var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
                Matrix2D shapeGradientsCartesian = 
                    jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
                Matrix2D deformation = BuildDeformationMatrix(shapeGradientsCartesian);

                // Contribution of this gauss point to the element stiffness matrix
                Matrix2D partial = deformation.Transpose() * (constitutive * deformation);
                double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight; //TODO: this is used by all methods that integrate. I should cache it.
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
            foreach (ElasticMaterial2D m in materialsAtGaussPoints) m.ClearState();
        }

        public void ClearMaterialStresses()
        {
            foreach (ElasticMaterial2D m in materialsAtGaussPoints) m.ClearStresses();
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
                foreach (ElasticMaterial2D material in materialsAtGaussPoints)
                    if (material.Modified) return true;
                return false;
            }
        }

        public void ResetMaterialModified()
        {
            foreach (ElasticMaterial2D material in materialsAtGaussPoints) material.ResetModified();
        }

        public void SaveMaterialState()
        {
            foreach (ElasticMaterial2D m in materialsAtGaussPoints) m.SaveState();
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
            IReadOnlyList<Matrix2D> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

            for (int gp = 0; gp < numGPs; ++gp)
            {
                IMatrix2D constitutive = materialsAtGaussPoints[gp].ConstitutiveMatrix;
                var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
                Matrix2D shapeGrandientsCartesian =
                    jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
                Matrix2D deformation = BuildDeformationMatrix(shapeGrandientsCartesian);

                strains[gp] = new double[3];
                deformation.Multiply(localDisplVector, strains[gp]);
                stresses[gp] = new double[3];
                constitutive.Multiply(new Vector(strains[gp]), stresses[gp]);
            }

            return (strains, stresses);
        }

        private Matrix2D BuildDeformationMatrix(Matrix2D shapeGradientsCartesian)
        {
            var deformation = new Matrix2D(3, 2 * Nodes.Count);
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                int col0 = 2 * nodeIdx;
                int col1 = 2 * nodeIdx + 1;

                deformation[0, col0] = shapeGradientsCartesian[nodeIdx, 0];
                deformation[1, col1] = shapeGradientsCartesian[nodeIdx, 1];
                deformation[2, col0] = shapeGradientsCartesian[nodeIdx, 1];
                deformation[2, col1] = shapeGradientsCartesian[nodeIdx, 0];
            }
            return deformation;
        }

        /// <summary>
        /// The shape function matrix is 2-by-2n, where n = is the number of shape functions. Row 0 corresponds to dof X, while
        /// row 1 to dof Y.
        /// </summary>
        private Matrix2D BuildShapeFunctionMatrix(Vector shapeFunctions)
        {
            var array2D = new double[2, 2 * shapeFunctions.Length];
            for (int i = 0; i < shapeFunctions.Length; ++i)
            {
                array2D[0, 2 * i] = shapeFunctions[i];
                array2D[1, 2 * i + 1] = shapeFunctions[i];
            }
            return new Matrix2D(array2D);
        }
    }
}
