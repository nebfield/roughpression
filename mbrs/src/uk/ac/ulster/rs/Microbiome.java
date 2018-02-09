package uk.ac.ulster.rs;

import rseslib.processing.classification.rules.roughset.RoughSetRuleClassifier;
import rseslib.structure.rule.Rule;
import rseslib.structure.rule.RuleWithStatistics;
import rseslib.structure.table.ArrayListDoubleDataTable;
import rseslib.structure.table.DoubleDataTable;
import rseslib.system.progress.StdOutProgress;

import java.io.*;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

public class Microbiome {

  public static void main(String[] args) throws Exception {

    if (args.length != 1) {
      System.out.println("Usage:");
      System.out.println("    java ... uk.ac.ulster.rs.Microbiome <data file>");
      System.out.println();
      System.out.println("Train a rough set rule classifier on microbiome census data");
      System.out.println("Outputs rules and rule statistics");
      System.exit(0);
    }

    DoubleDataTable table = new ArrayListDoubleDataTable(
        new File(args[0]),
        new StdOutProgress());

    // sanity check: verify the class labels
    System.out.println(table.attributes().nominalDecisionAttribute());
    RoughSetRuleClassifier rsc = new RoughSetRuleClassifier(null,
        table,
        new StdOutProgress());

    // output rules to text
    try(PrintWriter out = new PrintWriter("rules.txt")){
      out.println(rsc.getRules().toString());
    }

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
  }
}

