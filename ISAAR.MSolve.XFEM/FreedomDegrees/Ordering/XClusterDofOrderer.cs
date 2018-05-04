using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.FreedomDegrees.Ordering
{
    class XClusterDofOrderer
    {
        private readonly XCluster2D cluster;
        private readonly DofTable<DisplacementDof> constrainedDofs;
        private readonly DofTable<DisplacementDof> standardDofs;

        private XClusterDofOrderer(XCluster2D cluster, int numConstrainedDofs, DofTable<DisplacementDof> constrainedDofs, 
            int numStandardDofs, DofTable<DisplacementDof> standardDofs)
        {
            this.cluster = cluster;
            this.NumConstrainedDofs = numConstrainedDofs;
            this.constrainedDofs = constrainedDofs;
            this.NumStandardDofs = numStandardDofs;
            this.standardDofs = standardDofs;
        }

        public int NumConstrainedDofs { get; }
        public int NumStandardDofs { get; }

        public static XClusterDofOrderer CreateNodeMajor(Model2D model, XCluster2D cluster)
        {
            (int numConstrainedDofs, DofTable<DisplacementDof> constrainedDofs) = OrderConstrainedDofs(model.Constraints);
            (int numStandardDofs, DofTable<DisplacementDof> standardDofs) = OrderStandardDofs(model);

            int numTotalDofs = numStandardDofs;
            foreach (var subdomain in cluster.Subdomains)
            {
                subdomain.DofOrderer = XSubdomainDofOrderer.CreateNodeMajor(subdomain, numTotalDofs);
                numTotalDofs += subdomain.DofOrderer.NumEnrichedDofs;
            }
            return new XClusterDofOrderer(cluster, numConstrainedDofs, constrainedDofs, numStandardDofs, standardDofs);
        }

        public void ReorderStandardDofs(OrderingAMD orderingAlgorithm)
        {
            throw new NotImplementedException();
        }

        public void ReorderEnrichedSubdomainDofs(OrderingAMD orderingAlgorithm)
        {
            throw new NotImplementedException();
        }

        private static (int numConstrainedDofs, DofTable<DisplacementDof> constrainedDofs) OrderConstrainedDofs(
            ITable<XNode2D, DisplacementDof, double> constraints)
        {
            var constrainedDofs = new DofTable<DisplacementDof>();
            int counter = 0;
            foreach (Tuple<XNode2D, DisplacementDof, double> entry in constraints)
            {
                constrainedDofs[entry.Item1, entry.Item2] = counter++;
            }
            return (counter, constrainedDofs);
        }

        private static (int numStandardDofs, DofTable<DisplacementDof> standardDofs) OrderStandardDofs(Model2D model)
        {
            ITable<XNode2D, DisplacementDof, double> constraints = model.Constraints;
            var standardDofs = new DofTable<DisplacementDof>();
            int counter = 0;
            foreach (var node in model.Nodes)
            {
                if (!constraints.Contains(node, DisplacementDof.X)) standardDofs[node, DisplacementDof.X] = counter++;
                if (!constraints.Contains(node, DisplacementDof.Y)) standardDofs[node, DisplacementDof.Y] = counter++;
            }
            return (counter, standardDofs);
        }

        public void WriteToConsole()
        {
            Console.WriteLine("Standard free dofs for the whole domain:");
            Console.WriteLine(standardDofs);
            Console.WriteLine("Standard constrained dofs for the whole domain:");
            Console.WriteLine(constrainedDofs);
            for (int s = 0; s < cluster.Subdomains.Count; ++s)
            {
                Console.WriteLine($"Subdomain {s}:");
                cluster.Subdomains[s].DofOrderer.WriteToConsole();
            }
        }
    }
}
