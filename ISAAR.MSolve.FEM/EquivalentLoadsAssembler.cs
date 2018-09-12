using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM
{
    public static class EquivalentLoadsAssembler
    {
        public static IVector GetEquivalentNodalLoads(Subdomain subdomain, IElementMatrixProvider elementProvider, IVector solution, IVector dSolution) //TODOMaria this should also take as argument the nodal displacements of the constraints (after the refactoring)
        {
            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);
            var subdomainEquivalentNodalForces = new double[subdomain.Forces.Length];
            foreach (Element element in subdomain.ElementsDictionary.Values)
            {
                var isEmbeddedElement = element.ElementType is IEmbeddedElement;
                var elStart = DateTime.Now;
                IMatrix2D ElementK = elementProvider.Matrix(element);

                double[] localSolution = subdomain.CalculateElementNodalDisplacements(element, solution);
                double[] localdSolution = subdomain.CalculateElementNodalDisplacements(element, dSolution);

                var equivalentNodalForces = new double[localSolution.Length];
                ElementK.Multiply(new Vector(localSolution), equivalentNodalForces);
                subdomain.AddLocalVectorToGlobal(element, equivalentNodalForces, subdomainEquivalentNodalForces);

                times["addition"] += DateTime.Now - elStart;
            }

            var totalTime = DateTime.Now - totalStart;

            return new Vector(subdomainEquivalentNodalForces);
        }
    }
}
