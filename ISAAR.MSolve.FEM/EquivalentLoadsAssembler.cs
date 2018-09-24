using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM
{
    public class EquivalentLoadsAssembler : IEquivalentLoadsAssembler
    {
        private Subdomain subdomain;
        private Discretization.Interfaces.IElementMatrixProvider elementProvider;

        public EquivalentLoadsAssembler(Subdomain subdomain, Discretization.Interfaces.IElementMatrixProvider elementProvider)
        {
            this.subdomain = subdomain;
            this.elementProvider = elementProvider;
        }

        public IVector GetEquivalentNodalLoads(IVector solution, double constraintScalingFactor) 
        {
            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);
            var subdomainEquivalentNodalForces = new double[subdomain.Forces.Length];
            foreach (Element element in subdomain.ElementsDictionary.Values)
            {
                var isEmbeddedElement = element.ElementType is Interfaces.IEmbeddedElement;
                var elStart = DateTime.Now;
                IMatrix2D ElementK = elementProvider.Matrix(element);

                double[] localSolution = subdomain.CalculateElementNodalDisplacements(element, solution);
                double[] localdSolution = subdomain.CalculateElementIcrementalConstraintDisplacements(element, constraintScalingFactor);

                var equivalentNodalForces = new double[localSolution.Length];
                ElementK.Multiply(new Vector(localdSolution), equivalentNodalForces);
                subdomain.AddLocalVectorToGlobal(element, equivalentNodalForces, subdomainEquivalentNodalForces);

                times["addition"] += DateTime.Now - elStart;
            }

            var totalTime = DateTime.Now - totalStart;

            return new Vector(subdomainEquivalentNodalForces);
        }

        //private Vector CalculateElementPrescribedDisplacements(Element element)
        //{
        //    int numDofs = 0;
        //    IList<IList<DOFType>> nodalDofs = element.ElementType.GetElementDOFTypes(element);
        //    for (int nodeIdx = 0; nodeIdx < nodalDofs.Count; ++nodeIdx)
        //    {
        //        numDofs += nodalDofs[nodeIdx].Count;
        //    }

        //    //element.ElementType.GetElementDOFTypes(element) uses its own order for DOFTypes of the same node.
        //    //The vector u uses the same order as element.ElementType.GetElementDOFTypes(element).
        //    //node.Constraints may contain Constraint objects in different order of their DOFType
        //    //Perhaps use SortedSet instead of List for node.Constraints
        //    var u = new Vector(numDofs);
        //    for (int nodeIdx = 0; nodeIdx < nodalDofs.Count; ++nodeIdx)
        //    {
        //        Node node = element.Nodes[nodeIdx];
        //        if (node.Constraints.Count > 0)
        //        {
        //            foreach (var dof in nodalDofs[nodeIdx])
        //            {

        //            }
        //        }
        //    }
        //}
    }
}
