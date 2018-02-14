package uk.ac.ulster.rs;

import rseslib.processing.classification.rules.roughset.RoughSetRuleClassifier;
import rseslib.processing.discretization.DiscretizationFactory;
import rseslib.processing.transformation.TableTransformer;
import rseslib.processing.transformation.TransformationProvider;
import rseslib.processing.transformation.Transformer;
import rseslib.structure.rule.Rule;
import rseslib.structure.rule.RuleWithStatistics;
import rseslib.structure.table.ArrayListDoubleDataTable;
import rseslib.structure.table.DoubleDataTable;
import rseslib.system.Configuration;
import rseslib.system.progress.StdOutProgress;

import java.io.*;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;

public class Microbiome {

  public static void main(String[] args) throws Exception {

    if (args.length != 1) {
      System.out.println("Usage:");
      System.out.println("    java ... uk.ac.ulster.rs.Microbiome <data file>");
      System.out.println();
      System.out.println("Train a rough set rule classifier on microbiome census data");
      System.out.println("Outputs rules");
      System.exit(0);
    }

    DoubleDataTable table = new ArrayListDoubleDataTable(new File(args[0]),
        new StdOutProgress());

//    DoubleDataTable table = new ArrayListDoubleDataTable(new File("gut.arff"),
//        new StdOutProgress());

    // sanity check: verify the class labels
    System.out.println(table.attributes().nominalDecisionAttribute());

    // load default properties
    Properties rscProp = Configuration.loadDefaultProperties(RoughSetRuleClassifier.class);

    RoughSetRuleClassifier rsc = new RoughSetRuleClassifier(rscProp,
        table,
        new StdOutProgress());

    // store rules and support in hashmap
    Collection<Rule> rules = rsc.getRules();
    HashMap<String, Double> support = new HashMap<String, Double>();

    for (Rule r : rules) {
      support.put(r.toString(), ((RuleWithStatistics) r).getSupport());
    }

    // write hashmap to csv file
    Writer writer = new FileWriter("rule-support.csv");

    for (Map.Entry<String, Double> entry : support.entrySet()) {
      writer.append(entry.getKey())
        .append(",")
        .append(entry.getValue().toString());
    }
    writer.close();

    // output discretised data for analysis
    System.out.println("Discretisation used: " + rsc.getProperty("Discretization"));
    TransformationProvider discrProv = DiscretizationFactory.getDiscretizationProvider(rscProp);
    Transformer m_cDiscretizer = null;
    if (discrProv != null)
      m_cDiscretizer = discrProv.generateTransformer(table);
    DoubleDataTable discrTable = null;
    if (m_cDiscretizer != null)
      discrTable = TableTransformer.transform(table, m_cDiscretizer);
    if (discrTable != null)
      discrTable.storeArff(args[0], new File("discretised.arff"), new StdOutProgress());
  }

}

