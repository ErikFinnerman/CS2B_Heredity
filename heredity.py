import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():
    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):
                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """

    # Parse the input into cases for each person
    genes = dict()
    trait = dict()
    for person in people:
        if person in one_gene:
            genes[person] = 1
        elif person in two_genes:
            genes[person] = 2
        else:
            genes[person] = 0
        if person in have_trait:
            trait[person] = True
        else:
            trait[person] = False

    p_list = dict()
    for person in people:
        if people[person]['mother'] is None and people[person]['father'] is None: # Person does not have parents
            p_gene = PROBS['gene'][genes[person]]
            p_trait = PROBS['trait'][genes[person]][trait[person]]
            p = p_gene * p_trait
        else:  # If person has parents
            parent_genes = dict()
            p_parents = dict()  # Probability of getting one gene from the parent
            for parent in ('father', 'mother'):
                parent_genes[parent] = genes[people[person][parent]]
                if parent_genes[parent] == 1:
                    p_parents[parent] = 0.5 * (1 - PROBS['mutation'])+0.5*PROBS['mutation'] # 1 of the 2 genes is passed on randomly with prob = 0.5
                elif parent_genes[parent] == 2:
                    p_parents[parent] = 1 - PROBS['mutation']
                else: # I.e., parent does not have the gene
                    p_parents[parent] = PROBS['mutation']
            if genes[person] == 1:  # Gene either from mother or father, but not both
                p_gene = p_parents['father'] * (1 - p_parents['mother']) + p_parents['mother'] * (
                        1 - p_parents['father'])
            elif genes[person] == 2:  # Gene from both mother AND father
                p_gene = p_parents['father'] * p_parents['mother']
            elif genes[person] == 0:  # No gene from either mother OR father = NOT father AND NOT Mother
                p_gene = (1 - p_parents['father']) * (1 - p_parents['mother'])

            p_trait = PROBS['trait'][genes[person]][trait[person]]
            p = p_gene * p_trait
        p_list[person] = p
    joint_p = 1
    for p in p_list.values():
        joint_p = joint_p * p
    return joint_p


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """

    # Update for each person
    for person in probabilities:

        #Add probabilities for one gene, 2 genes or no genes
        if person in one_gene:
            probabilities[person]["gene"][1] += p
        elif person in two_genes:
            probabilities[person]["gene"][2] += p
        else:
            probabilities[person]["gene"][0] += p

        # Add probabilities for Trait
        if person in have_trait:
            probabilities[person]['trait'][True] += p
        else:
            probabilities[person]['trait'][False] += p

def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """

    #Make the update for each person
    for person in probabilities:

    # Update 'gene' probability distribution
        # Get the sum of all gene values
        gene_sum=sum(probabilities[person]['gene'].values())
        if gene_sum==0:
            continue
        else:
            probabilities[person]['gene'][1]=probabilities[person]['gene'][1]/gene_sum
            probabilities[person]['gene'][2] = probabilities[person]['gene'][2] / gene_sum
            probabilities[person]['gene'][0] = probabilities[person]['gene'][0] / gene_sum

    # Update 'trait' probability distribution
        # Get the sum of all gene values
        trait_sum = sum(probabilities[person]['trait'].values())
        if trait_sum==0:
            continue
        else:
            probabilities[person]['trait'][True] = probabilities[person]['trait'][True] / gene_sum
            probabilities[person]['trait'][False] = probabilities[person]['trait'][False] / gene_sum

if __name__ == "__main__":
    main()
    test_people={'Harry': {'name': 'Harry', 'mother': 'Lily', 'father': 'James', 'trait': None},
              'James': {'name': 'James', 'mother': None, 'father': None, 'trait': True},
              'Lily': {'name': 'Lily', 'mother': None, 'father': None, 'trait': False}}
    joint_probability(test_people, {"Harry"}, {"James"}, {"James"})


