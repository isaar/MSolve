using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Logging
{
    public enum NEUOutputVector
    {
        TranslationX = 2,
        TranslationY,
        TranslationZ,
        RotationX = 6,
        RotationY,
        RotationZ
    }

    public class NEUOutputVectorProperties
    {
        private readonly int id;
        private readonly string title;
        private readonly string componentVector;
        private readonly IDofType dofType;

        public int ID { get { return id; } }
        public string Title { get { return title; } }
        public string ComponentVector { get { return componentVector; } }
        public IDofType DOFType { get { return dofType; } }

        public NEUOutputVectorProperties(int id, string title, string componentVector, IDofType dofType)
        {
            this.id = id;
            this.title = title;
            this.componentVector = componentVector;
            this.dofType = dofType;
        }
    }

    public class NEUWriter
    {
        private const int outputSet = 1;
        private const string BlockMarker = "   -1";
        private readonly Dictionary<NEUOutputVector, NEUOutputVectorProperties> vectorTitles = new Dictionary<NEUOutputVector, NEUOutputVectorProperties>()
        {
            { NEUOutputVector.TranslationX, new NEUOutputVectorProperties(2, "X Translation", "2,0,0", StructuralDof.TranslationX) },
            { NEUOutputVector.TranslationY, new NEUOutputVectorProperties(3, "Y Translation", "0,3,0", StructuralDof.TranslationY) },
            { NEUOutputVector.TranslationZ, new NEUOutputVectorProperties(4, "Z Translation", "0,0,4", StructuralDof.TranslationZ) },
            { NEUOutputVector.RotationX, new NEUOutputVectorProperties(6, "X Rotation", "6,0,0", StructuralDof.RotationX) },
            { NEUOutputVector.RotationY, new NEUOutputVectorProperties(7, "Y Rotation", "0,7,0", StructuralDof.RotationY) },
            { NEUOutputVector.RotationZ, new NEUOutputVectorProperties(8, "Z Rotation", "0,0,8", StructuralDof.RotationZ) }
        };

        private readonly ILinearSystem subdomain;
        private readonly Model model;

        public NEUWriter(Model model, ILinearSystem subdomain)
        {
            this.model = model;
            this.subdomain = subdomain;
        }

        private List<string> MakeOutputSet(string title)
        {
            var lines = new List<string>();
            lines.Add(BlockMarker);
            lines.Add("   450");
            lines.Add(String.Format("{0},", outputSet));
            lines.Add(title);
            lines.Add(String.Format("{0},{1},", 0, 1));
            lines.Add("0.0,");
            lines.Add("1,");
            lines.Add("<No comments>");
            lines.Add(BlockMarker);

            return lines;
        }

        private List<string> MakeNodalDataVector(NEUOutputVectorProperties vectorProperties)
        {
            var lines = new List<string>();
            lines.Add(String.Format("{0},{1},1,", outputSet, vectorProperties.ID));
            lines.Add(vectorProperties.Title);
            lines.Add(String.Format("{0},{1},{2},", 0, 0, 0));
            lines.Add(String.Format("{0},0,0,0,0,0,0,0,", vectorProperties.ComponentVector));
            lines.Add("0,0,0,0,0,0,0,0,0,0,");
            lines.Add("0,");
            lines.Add(String.Format("{0},{1},1,7,", 0, 0));
            lines.Add("0,1,1,");

            var embeddedNodeValues = CalculateEmbeddedNodeValues();
            foreach (var node in model.Nodes)
            {
                var key = new Tuple<int, IDofType>(node.ID, vectorProperties.DOFType);
                if (embeddedNodeValues.ContainsKey(key))
                    lines.Add(String.Format(CultureInfo.InvariantCulture, "{0},{1},", node.ID, embeddedNodeValues[key]));
                else
                {
                    bool nodeExists = model.GlobalDofOrdering.GlobalFreeDofs.TryGetDataOfRow(node,
                        out IReadOnlyDictionary<IDofType, int> dofTypesIndices);
                    if (nodeExists)
                    {
                        if (!dofTypesIndices.ContainsKey(vectorProperties.DOFType)) continue;
                        lines.Add(String.Format(CultureInfo.InvariantCulture, "{0},{1},", 
                            node.ID, subdomain.Solution[dofTypesIndices[vectorProperties.DOFType]]));
                    }
                        
                }
               
            }
            lines.Add("-1,0.,");

            return lines;
        }

        private Dictionary<Tuple<int, IDofType>, double> CalculateEmbeddedNodeValues()
        {
            var embeddedNodeValues = new Dictionary<Tuple<int, IDofType>, double>();

            foreach (var element in model.Elements)
            {
                IEmbeddedElement e = element.ElementType as IEmbeddedElement;
                if (e == null) continue;

                var superElementNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
                var superElementDOFs = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
                var superElementVector = new double[superElementDOFs.SelectMany(x => x).Count()];

                int index = 0;
                for (int i = 0; i < superElementDOFs.Count; i++)
                    for (int j = 0; j < superElementDOFs[i].Count; j++)
                    {
                        bool isDofFree = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(
                            superElementNodes[i], superElementDOFs[i][j], out int dof);
                        if (isDofFree) superElementVector[index] = subdomain.Solution[dof];
                        index++;
                    }

                index = 0;
                var elementVector = e.GetLocalDOFValues(element, superElementVector);
                var elementDOFs = element.ElementType.GetElementDofTypes(element);
                for (int i = 0; i < elementDOFs.Count; i++)
                    for (int j = 0; j < elementDOFs[i].Count; j++)
                    {
                        var key = new Tuple<int, IDofType>(element.Nodes[i].ID, elementDOFs[i][j]);
                        //int dof = model.NodalDOFsDictionary[element.Nodes[i].ID][elementDOFs[i][j]];

                        if (!embeddedNodeValues.ContainsKey(key))
                            embeddedNodeValues.Add(key, elementVector[index]);
                        index++;
                    }

                //var embeddedDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypesForDOFEnumeration(element);
                //foreach (var embeddedNode in e.EmbeddedNodes)
                //{
                //    Element hostElement = embeddedNode.EmbeddedInElement;
                //    var hostDOFTypes = hostElement.ElementType.DOFEnumerator.GetDOFTypes(hostElement);
                //    var hostValues = new double[hostDOFTypes.SelectMany(x => x).Count()];
                //    int index = 0;
                //    for (int i = 0; i < hostDOFTypes.Count; i++)
                //        for (int j = 0; j < hostDOFTypes[i].Count; j++)
                //        {
                //            int dof = model.NodalDOFsDictionary[hostElement.Nodes[i].ID][hostDOFTypes[i][j]];
                //            if (dof > -1)
                //                hostValues[index] = subdomain.Solution[dof];
                //            index++;
                //        }

                //    var embeddedValues = e.GetLocalDOFValues(hostElement, hostValues);
                //    int nodeIndex = element.Nodes.IndexOf(embeddedNode.Node);
                //    index = 0;
                //    for (int i = 0; i < nodeIndex; i++)
                //        index += embeddedDOFTypes[nodeIndex].Count;
                //    //for (int i = 0; i < embeddedDOFTypes.Count; i++)
                //    for (int j = 0; j < embeddedDOFTypes[nodeIndex].Count; j++)
                //        {
                //            var key = new Tuple<int, DOFType>(element.Nodes[nodeIndex].ID, embeddedDOFTypes[nodeIndex][j]);
                //            int dof = model.NodalDOFsDictionary[element.Nodes[nodeIndex].ID][embeddedDOFTypes[nodeIndex][j]];

                //            if (dof < 0 && !embeddedNodeValues.ContainsKey(key))
                //                embeddedNodeValues.Add(key, embeddedValues[index]);
                //            index++;
                //        }
                //}
            }

            return embeddedNodeValues;
        }

        private List<string> MakeDummyTotalNodalVectors()
        {
            var lines = new List<string>();
            lines.Add(String.Format("{0},{1},1,", outputSet, 1));
            lines.Add("Total Translation");
            lines.Add(String.Format("{0},{1},{2},", 0, 0, 0));
            lines.Add("2,3,4,0,0,0,0,0,0,0,");
            lines.Add("0,0,0,0,0,0,0,0,0,0,");
            lines.Add("0,");
            lines.Add(String.Format("{0},{1},1,7,", 0, 0));
            lines.Add("1,1,1,");
            foreach (var node in model.Nodes)
            {
                double translation = 0;
                bool nodeExists = model.GlobalDofOrdering.GlobalFreeDofs.TryGetDataOfRow(node, 
                    out IReadOnlyDictionary<IDofType, int> dofTypesIndices);
                if (nodeExists)
                {
                    foreach (var dof in dofTypesIndices.Where(
                        x => x.Key == StructuralDof.TranslationX || x.Key == StructuralDof.TranslationY || x.Key == StructuralDof.TranslationZ).Select(x => x.Value))
                    {
                        translation += Math.Pow(subdomain.Solution[dof], 2);
                    }
                }
                lines.Add(String.Format(CultureInfo.InvariantCulture, "{0},{1},", node.ID, Math.Sqrt(translation)));
            }
            lines.Add("-1,0.,");

            lines.Add(String.Format("{0},{1},1,", outputSet, 5));
            lines.Add("Total Rotation");
            lines.Add(String.Format("{0},{1},{2},", 0, 0, 0));
            lines.Add("6,7,8,0,0,0,0,0,0,0,");
            lines.Add("0,0,0,0,0,0,0,0,0,0,");
            lines.Add("0,");
            lines.Add(String.Format("{0},{1},1,7,", 0, 0));
            lines.Add("1,1,1,");
            foreach (var node in model.Nodes)
            {
                double translation = 0;
                bool nodeExists = model.GlobalDofOrdering.GlobalFreeDofs.TryGetDataOfRow(node, 
                    out IReadOnlyDictionary<IDofType, int> dofTypesIndices);
                if (nodeExists)
                {
                    foreach (var dof in dofTypesIndices.Where(
                        x => x.Key == StructuralDof.RotationX || x.Key == StructuralDof.RotationY || x.Key == StructuralDof.RotationZ).Select(x => x.Value))
                    {
                        translation += Math.Pow(subdomain.Solution[dof], 2);
                    }
                }
                lines.Add(String.Format(CultureInfo.InvariantCulture, "{0},{1},", node.ID, Math.Sqrt(translation)));
            }
            lines.Add("-1,0.,");

            return lines;
        }

        private List<string> MakeNodalDataVectors()
        {
            var lines = new List<string>();
            lines.Add(BlockMarker);
            lines.Add("  1051");
            lines.AddRange(MakeDummyTotalNodalVectors());
            foreach (var prop in vectorTitles)
                lines.AddRange(MakeNodalDataVector(prop.Value));
            lines.Add(BlockMarker);

            return lines;
        }

        private List<string> MakeHeader()
        {
            var lines = new List<string>();
            lines.Add(BlockMarker);
            lines.Add("   100");
            lines.Add("<NULL>");
            lines.Add("10.2,");
            lines.Add(BlockMarker);
            return lines;
        }

        public void WriteToFile(string fileName)
        {
            var lines = new List<string>();
            lines.AddRange(MakeHeader());
            lines.AddRange(MakeOutputSet("Output set 1"));
            lines.AddRange(MakeNodalDataVectors());

            File.WriteAllLines(fileName, lines);
        }
    }
}
