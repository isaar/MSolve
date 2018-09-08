using ISAAR.MSolve.XFEM.FreedomDegrees;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.XFEM.Entities
{
    class ContinuityEquations
    {
        public ContinuityEquations(int numEquations, IReadOnlyList<XNode2D> nodes, IReadOnlyList<EnrichedDof> dofs,
            IReadOnlyList<XSubdomain2D> positiveSubdomains, IReadOnlyList<XSubdomain2D> negativeSubdomains)
        {
            this.NumEquations = numEquations;
            this.Nodes = nodes;
            this.Dofs = dofs;
            this.PositiveSubdomains = positiveSubdomains;
            this.NegativeSubdomains = negativeSubdomains;
        }

        public enum EquationOrder { DofMajor, SubdomainMajor};

        public int NumEquations { get; }
        public IReadOnlyList<XNode2D> Nodes { get; }
        public IReadOnlyList<EnrichedDof> Dofs { get; }
        public IReadOnlyList<XSubdomain2D> PositiveSubdomains { get; }
        public IReadOnlyList<XSubdomain2D> NegativeSubdomains { get; }


        public class Builder
        {
            private readonly Dictionary<XNode2D, Dictionary<EnrichedDof, SortedSet<XSubdomain2D>>> data;

            public Builder()
            {
                this.data = new Dictionary<XNode2D, Dictionary<EnrichedDof, SortedSet<XSubdomain2D>>>();
            }

            public ContinuityEquations Build(EquationOrder order)
            {
                var nodes = new List<XNode2D>();
                var dofs = new List<EnrichedDof>();
                var positiveSubdomains = new List<XSubdomain2D>();
                var negativeSubdomains = new List<XSubdomain2D>();

                if (order == EquationOrder.DofMajor)
                {
                    foreach (var nodeData in data)
                    {
                        XNode2D node = nodeData.Key;
                        foreach (var dofData in nodeData.Value)
                        {
                            EnrichedDof dof = dofData.Key;
                            if (dofData.Value.Count > 1) // There must be >= 2 subdomains, in order to enforce continuity
                            {
                                XSubdomain2D[] subdomains = dofData.Value.ToArray();
                                for (int i = 0; i < subdomains.Length - 1; ++i)
                                {
                                    nodes.Add(node);
                                    dofs.Add(dof);
                                    positiveSubdomains.Add(subdomains[i]);
                                    negativeSubdomains.Add(subdomains[i + 1]);
                                }
                            }
                        }
                    }
                    return new ContinuityEquations(nodes.Count, nodes, dofs, positiveSubdomains, negativeSubdomains);
                }
                else
                {
                    throw new NotImplementedException();
                }
            }

            public void Register(XNode2D node, EnrichedDof dof, XSubdomain2D subdomain)
            {
                bool nodeExists = data.TryGetValue(node, out Dictionary<EnrichedDof, SortedSet<XSubdomain2D>> nodeData);
                if (!nodeExists)
                {
                    nodeData = new Dictionary<EnrichedDof, SortedSet<XSubdomain2D>>();
                    data.Add(node, nodeData);
                }

                bool dofExists = nodeData.TryGetValue(dof, out SortedSet<XSubdomain2D> subdomains);
                if (!dofExists)
                {
                    subdomains = new SortedSet<XSubdomain2D>();
                    nodeData.Add(dof, subdomains);
                }

                subdomains.Add(subdomain);
            }

        }
    }
}
