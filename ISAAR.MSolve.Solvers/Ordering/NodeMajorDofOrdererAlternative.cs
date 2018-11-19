using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Free dofs are assigned global (actually subdomain) indices in a node major fashion: The dofs of the first node are 
    /// numbered, then the dofs of the second node, etc. Constrained dofs are ignored.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class NodeMajorDofOrdererAlternative: IDofOrderer
    {
        public IDofOrdering OrderDofs(IStructuralModel_v2 model, ISubdomain_v2 subdomain)
        {
            int totalDOFs = 0;
            Dictionary<int, List<DOFType>> nodalDOFTypesDictionary = new Dictionary<int, List<DOFType>>(); //TODO: use Set isntead of List
            foreach (IElement element in subdomain.ΙElementsDictionary.Values)
            {
                for (int i = 0; i < element.INodes.Count; i++)
                {
                    if (!nodalDOFTypesDictionary.ContainsKey(element.INodes[i].ID))
                        nodalDOFTypesDictionary.Add(element.INodes[i].ID, new List<DOFType>());
                    nodalDOFTypesDictionary[element.INodes[i].ID].AddRange(element.IElementType.DOFEnumerator.GetDOFTypesForDOFEnumeration(element)[i]);
                }
            }

            var freeDofs = new DofTable();
            foreach (INode node in subdomain.Nodes)
            {
                //List<DOFType> dofTypes = new List<DOFType>();
                //foreach (Element element in node.ElementsDictionary.Values)
                //{
                //    if (elementsDictionary.ContainsKey(element.ID))
                //    {
                //        foreach (DOFType dof in element.ElementType.DOFTypes)
                //            dofTypes.Add(dof);
                //    }
                //}

                Dictionary<DOFType, int> dofsDictionary = new Dictionary<DOFType, int>();
                //foreach (DOFType dofType in dofTypes.Distinct<DOFType>())
                foreach (DOFType dofType in nodalDOFTypesDictionary[node.ID].Distinct<DOFType>())
                {
                    int dofID = 0;
                    #region removeMaria
                    //foreach (DOFType constraint in node.Constraints)
                    //{
                    //    if (constraint == dofType)
                    //    {
                    //        dofID = -1;
                    //        break;
                    //    }
                    //}
                    #endregion

                    foreach (var constraint in ((FEM.Entities.Node)node).Constraints) //TODO: access the constraints from the subdomain
                    {
                        if (constraint.DOF == dofType)
                        {
                            dofID = -1;
                            break;
                        }
                    }

                    //var embeddedNode = embeddedNodes.Where(x => x.Node == node).FirstOrDefault();
                    ////if (node.EmbeddedInElement != null && node.EmbeddedInElement.ElementType.GetDOFTypes(null)
                    ////    .SelectMany(d => d).Count(d => d == dofType) > 0)
                    ////    dofID = -1;
                    //if (embeddedNode != null && embeddedNode.EmbeddedInElement.ElementType.DOFEnumerator.GetDOFTypes(null)
                    //    .SelectMany(d => d).Count(d => d == dofType) > 0)
                    //    dofID = -1;

                    if (dofID == 0)
                    {
                        dofID = totalDOFs;
                        freeDofs[node, dofType] = dofID;
                        totalDOFs++;
                    }
                }
            }

            return new FreeDofOrdering(totalDOFs, freeDofs, model.NodalDOFsDictionary);
        }
    }
}
