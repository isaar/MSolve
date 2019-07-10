using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

//TODO: decide on consistent collections, indexes/keys and naming
namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    public interface ICrackDescription
    {
        /// <summary>
        /// All nodes enriched with Heaviside.
        /// </summary>
        IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesAll { get; }

        /// <summary>
        /// Nodes that were first enriched with Heaviside in the current iteration. Subset of <see cref="crackBodyNodesAll"/>
        /// </summary>
        IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesNew { get; }

        /// <summary>
        /// Nodes that were already enriched with Heaviside, but their crack body level set changes in the current iteration.
        /// Most of the time this should happen to nodes of the element containing the crack tip, thus this set should be empty.
        /// </summary>
        IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesModified { get; }

        /// <summary>
        /// Nodes that were already enriched with Heaviside, but belong to elements, where at least one node belongs to one of
        /// the sets <see cref="crackBodyNodesNew"/>, <see cref="crackBodyNodesModified"/>, <see cref="crackTipNodesNew"/> or 
        /// <see cref="crackTipNodesOld"/>. The stiffness matrix of these elements will change, causing a change in the stiffness
        /// terms associated with these nodes.
        /// </summary>
        IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesNearModified { get;}

        /// <summary>
        /// Nodes that belong to elements intersected by the 0 crack body level set, but must not be enriched with Heaviside to
        /// avoid singularity in the global stiffness matrix.
        /// </summary>
        IReadOnlyDictionary<CrackBodyEnrichment2D, ISet<XNode>> CrackBodyNodesRejected { get; }

        IReadOnlyList<CartesianPoint> CrackTips { get; }
        IReadOnlyDictionary<CartesianPoint, IReadOnlyList<XContinuumElement2D>> CrackTipElements { get; }
        IReadOnlyDictionary<CartesianPoint, IPropagator> CrackTipPropagators { get; }

        /// <summary>
        /// Nodes that are currently enriched with tip functions.
        /// </summary>
        IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode>> CrackTipNodesNew { get; }

        /// <summary>
        /// Nodes that were enriched with tip functions in the previous iteration, but now are not. 
        /// </summary>
        IReadOnlyDictionary<CrackTipEnrichments2D, ISet<XNode>> CrackTipNodesOld { get; }

        /// <summary>
        /// Elements with at least one node whose enrichment has changed (added Heaviside, tip or removed tip functions) since 
        /// the previous iteration. These are the only elements whose stiffness matrices are modified.
        /// </summary>
        ISet<XContinuumElement2D> ElementsModified { get; }

        IReadOnlyList<IEnrichmentItem2D> Enrichments { get; }
        BidirectionalMesh2D<XNode, XContinuumElement2D> Mesh { get; } //TODO: abstract this

        IReadOnlyList<ISingleCrack> SingleCracks { get; } //Not sure about this one
        
        void Propagate(Dictionary<int, Vector> totalFreeDisplacements);
        void UpdateEnrichments();
    }
}
