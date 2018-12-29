using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Preprocessor.Meshes;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;

namespace ISAAR.MSolve.Preprocessor.UI
{
    /// <summary>
    /// Utility class to facilitate the definition of the model: nodes, elements, boundary conditions, loads, etc.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PreprocessorModel
    {
        internal enum ProblemDimensions //TODO: what about beams, trusses, plates and shells?
        {
            TwoDimensionalPlaneStress, TwoDimensionalPlaneStrain, ThreeDimensional
        }

        private readonly ElasticMaterial2D planeStrainMaterial;
        private readonly double thickness;
        private readonly FEM.Entities.Model model;
        private bool isAssembled;

        private PreprocessorModel(ProblemDimensions dimensions, double thickness, ElasticMaterial2D planeStrainMaterial)
        {
            this.isAssembled = false;
            this.planeStrainMaterial = planeStrainMaterial;
            this.thickness = thickness;
            this.model = new FEM.Entities.Model();
            this.model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 }); //TODO: let the user decide how many subdomains there will be
            this.Dimensions = dimensions;
        }

        /// <summary>
        /// Initializes a model for 2D plane stress problems.
        /// </summary>
        /// <param name="thickness">The 3rd dimension of the domain that is assumed to be much smaller than the main 2.</param>
        /// <returns></returns>
        public static PreprocessorModel Create2DPlaneStress(double thickness)
        {
            return new PreprocessorModel(ProblemDimensions.TwoDimensionalPlaneStress, thickness, null);
        }

        /// <summary>
        /// Initializes a model for 2D plane strain problems.
        /// </summary>
        /// <param name="planeStrainMaterial">The material properties for a plane strain problem. For now they must be the same  
        ///     for the whole domain.</param>
        /// <returns></returns>
        public static PreprocessorModel Create2DPlaneStrain(ElasticMaterial2D planeStrainMaterial) //TODO: extend this for non linear materials
        {
            if (planeStrainMaterial.StressState != StressState2D.PlaneStrain)
            {
                throw new ArgumentException("The provided material is not for plane strain problems");
            }
            return new PreprocessorModel(ProblemDimensions.TwoDimensionalPlaneStrain, 1.0, planeStrainMaterial);
        }

        /// <summary>
        /// Initializes a model for 3D problems. Not fully supported yet.
        /// </summary>
        /// <returns></returns>
        public static PreprocessorModel Create3D()
        {
            return new PreprocessorModel(ProblemDimensions.ThreeDimensional, 0.0, null);
        }

        internal Model CoreModel
        {
            get
            {
                if (!isAssembled)
                {
                    model.ConnectDataStructures();
                    isAssembled = true;
                }
                return model;
            }
        }

        internal ProblemDimensions Dimensions { get; }

        internal ElasticMaterial2D PlainStrainMaterial
        {
            get
            {
                if (Dimensions == ProblemDimensions.TwoDimensionalPlaneStrain) return planeStrainMaterial;
                else throw new InvalidOperationException("This property can be used only for 2D plane strain problems");
            }
        }

        internal double TimeStep { get; private set; }

        internal double TotalDuration { get; private set; }

        /// <summary>
        /// Adds a new node to the model.
        /// </summary>
        /// <param name="node">The node's coordinates and ID.</param>
        public void AddNode(Node node)
        {
            int numNodes = model.NodesDictionary.Count;
            model.NodesDictionary.Add(numNodes, node);
        }

        /// <summary>
        /// Adds a new finite element to the model.
        /// </summary>
        /// <param name="element">The element type (e.g. continuum element, beam, etc.)</param>
        /// <param name="elementNodes">The element's nodes. Beware of their order.</param>
        public void AddElement(IFiniteElement element, IReadOnlyList<Node> elementNodes)
        {
            int numElementsTotal = model.Elements.Count;
            int numElementsSubdomain = model.Subdomains[0].ElementsDictionary.Count;

            var elementWrapper = new Element() { ID = numElementsTotal, ElementType = element };
            foreach (Node node in elementNodes) elementWrapper.AddNode(node);
            model.ElementsDictionary.Add(numElementsTotal, elementWrapper);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(numElementsSubdomain, elementWrapper); //TODO: let the user decide which subdomain it will be added to.
        }

        /// <summary>
        /// Adds a whole mesh, which is defined by its nodes and elements, to the model.
        /// </summary>
        /// <param name="nodes">The nodes of the mesh.</param>
        /// <param name="elements">The elements of the mesh.</param>
        /// <param name="material">The common material of all elements in the mesh.</param>
        /// <param name="dynamicProperties">Optional material properties for dynamic analysis. Common for all elements in the 
        ///     mesh.</param>
        public void AddMesh2D(IReadOnlyList<Node2D> nodes, IReadOnlyList<CellConnectivity2D> elements, 
            IFiniteElementMaterial material, DynamicMaterial dynamicProperties = null)
        {
            if (Dimensions == ProblemDimensions.ThreeDimensional)
            {
                throw new InvalidOperationException(
                    "This method can be used only for 2D problems (plane strain or plane stress)");
            }

            // Nodes
            int numNodesCurrent = model.NodesDictionary.Count;
            for (int i = 0; i < nodes.Count; ++i) model.NodesDictionary.Add(numNodesCurrent + i, nodes[i]);

            // Elements
            int numElementsCurrent = model.Elements.Count;
            int numElementsSubdomain = model.Subdomains[0].ElementsDictionary.Count;
            var factory = new ContinuumElement2DFactory(thickness, (ElasticMaterial2D)material, dynamicProperties); //TODO: extend the factory to other materials
            for (int i = 0; i < elements.Count; ++i)
            {
                ContinuumElement2D element = factory.CreateElement(elements[i].CellType, elements[i].Vertices);
                var elementWrapper = new Element() { ID = numElementsCurrent + i, ElementType = element };
                foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(numElementsCurrent + i, elementWrapper);
                model.SubdomainsDictionary[0].ElementsDictionary.Add(numElementsSubdomain + i, elementWrapper); //TODO: let the user decide which subdomain it will be added to.
            }
        }

        /// <summary>
        /// Enforces a displacement for the provided node.
        /// </summary>
        /// <param name="node"></param>
        /// <param name="freedomDegree">The axis of the prescribed displacement.</param>
        /// <param name="signedMagnitude">The magnitude of the prescribed displacement with its sign.</param>
        public void ApplyPrescribedDisplacement(Node node, DOFType freedomDegree, double signedMagnitude)
        {
            if (signedMagnitude == 0.0)
            {
                node.Constraints.Add(new Constraint() { DOF = freedomDegree, Amount = 0.0 });
            }
            else throw new NotImplementedException("For now only 0 displacements are supported");
        }

        /// <summary>
        /// Enforces a displacement for all provided nodes.
        /// </summary>
        /// <param name="nodes"></param>
        /// <param name="freedomDegree">The axis of the prescribed displacement.</param>
        /// <param name="signedMagnitude">The magnitude of the prescribed displacement with its sign.</param>
        public void ApplyPrescribedDisplacements(IEnumerable<Node> nodes, DOFType freedomDegree, double signedMagnitude)
        {
            foreach (Node node in nodes) ApplyPrescribedDisplacement(node, freedomDegree, signedMagnitude);
        }

        /// <summary>
        /// Enforces a displacement for all provided nodes.
        /// </summary>
        /// <param name="nodeSelector">A method that returns true for the nodes to be selected and false for the rest.</param>
        /// <param name="freedomDegree">The axis of the prescribed displacement.</param>
        /// <param name="signedMagnitude">The magnitude of the prescribed displacement with its sign.</param>
        public void ApplyPrescribedDisplacements(Func<Node, bool> nodeSelector, DOFType freedomDegree, double signedMagnitude)
        {
            ApplyPrescribedDisplacements(model.Nodes.Where(nodeSelector), freedomDegree, signedMagnitude);
        }

        /// <summary>
        /// Applies a load at the provided node.
        /// </summary>
        /// <param name="node"></param>
        /// <param name="freedomDegree">The axis of the nodal load.</param>
        /// <param name="signedMagnitude">The magnitude of the nodal load with its sign.</param>
        public void ApplyNodalLoad(Node node, DOFType freedomDegree, double signedMagnitude)
        {
            model.Loads.Add(new Load() { Amount = signedMagnitude, Node = node, DOF = freedomDegree });
        }

        /// <summary>
        /// Applies the provided load to all provided nodes.
        /// </summary>
        /// <param name="nodes"></param>
        /// <param name="freedomDegree">The axis of the nodal load.</param>
        /// <param name="signedMagnitudeIndividual">The magnitude of the nodal load with its sign.</param>
        public void ApplyNodalLoads(IEnumerable<Node> nodes, DOFType freedomDegree, double signedMagnitudeIndividual)
        {
            foreach (Node node in nodes) ApplyNodalLoad(node, freedomDegree, signedMagnitudeIndividual);
        }

        /// <summary>
        /// Applies the provided load to all provided nodes.
        /// </summary>
        /// <param name="nodes">A method that returns true for the nodes to be selected and false for the rest.</param>
        /// <param name="freedomDegree">The axis of the nodal load.</param>
        /// <param name="signedMagnitudeIndividual">The magnitude of the nodal load with its sign.</param>
        public void ApplyNodalLoads(Func<Node, bool> nodeSelector, DOFType freedomDegree, double signedMagnitudeIndividual)
        {
            ApplyNodalLoads(model.Nodes.Where(nodeSelector), freedomDegree, signedMagnitudeIndividual);
        }

        /// <summary>
        /// Divides the provided load by the number of provided nodes and then applies it to all of them.
        /// </summary>
        /// <param name="nodes"></param>
        /// <param name="freedomDegree">The axis of the nodal load.</param>
        /// <param name="signedMagnitudeTotal">The magnitude of the nodal load with its sign.</param>
        public void DistributeNodalLoadEvenly(IEnumerable<Node> nodes, DOFType freedomDegree, double signedMagnitudeTotal)
        {
            int numNodes = nodes.Count();
            double individualLoad = signedMagnitudeTotal / numNodes;
            foreach (Node node in nodes) ApplyNodalLoad(node, freedomDegree, individualLoad);
        }

        /// <summary>
        /// Divides the provided load by the number of provided nodes and then applies it to all of them.
        /// </summary>
        /// <param name="nodes">A method that returns true for the nodes to be selected and false for the rest.</param>
        /// <param name="freedomDegree">The axis of the nodal load.</param>
        /// <param name="signedMagnitudeTotal">The magnitude of the nodal load with its sign.</param>
        public void DistributeNodalLoadEvenly(Func<Node, bool> nodeSelector, DOFType freedomDegree, double signedMagnitudeTotal)
        {
            DistributeNodalLoadEvenly(model.Nodes.Where(nodeSelector), freedomDegree, signedMagnitudeTotal);
        }

        /// <summary>
        /// Applies a ground motion to the whole structure. Only 1 can be active, but it may be applied to various freedom 
        /// degrees with various magnifications.
        /// </summary>
        /// <param name="accelerogramPath">The absolute path of the file that contains the ground acceleration at each 
        ///     timestep.</param>
        /// <param name="magnificationFactors">The acceleration at along each freedom degree will be equal to the value read 
        ///     from <paramref name="accelerogramPath"/> times the magnification factor for that degree. Do not include freedom 
        ///     degrees with 0 magnifications.</param>
        /// <param name="timeStep">The time between two successive accelerations read from 
        ///     <paramref name="accelerogramPath"/>.</param>
        /// <param name="TotalDuration">The total duration of the ground motion.</param>
        /// <param name="magnificationFactor"></param>
        public void SetGroundMotion(string accelerogramPath, Dictionary<DOFType, double> magnificationFactors, double timeStep,
            double TotalDuration)
        {
            this.TimeStep = timeStep;
            this.TotalDuration = TotalDuration;

            model.MassAccelerationHistoryLoads.Clear();
            foreach (var dofMagnificationPair in magnificationFactors)
            {
                model.MassAccelerationHistoryLoads.Add(
                    new MassAccelerationHistoryLoad(accelerogramPath, dofMagnificationPair.Value)
                    {
                        DOF = dofMagnificationPair.Key
                    });
            }
        }
    }
}
