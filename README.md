# Ethereum Foundation Initiatives

Elliptic curve algorithms (ECDSA, Ed25519) are vulnerable to Shor’s algorithm that runs on a quantum computer. According to estimates, such attacks may become feasible under 12 minutes by 2030.

Most people will migrate their ECDSA/Ed25519 accounts to post-quantum accounts. But some will be inevitably left over. For example, Estonian banker Rain Lohmus lost the private key to his account, containing 250 000 ETH. A quantum adversary could steal his money by deriving the account private key from the public key. This is not a single case — there are such examples for every major blockchain (e.g. Satoshi’s possibly dormant Bitcoin accounts).

Most chains would likely “freeze” legacy ECDSA/Ed25519 accounts to prevent theft, devaluation and loss of reputation. But in order to avoid censorship, chains will have to come up with ways for users to recover accounts after the freeze. The most likely way would perhaps be to generate a zero-knowledge proof, stating “I know a secret phrase S, such that my account with address A can be derived from it.”

Additionally, the consensus-level algorithm BLS is also vulnerable. It’s replacements XLS and XMSS have huge signatures. Additionally, BLS can aggregate thousands of validator signatures into one, whereas XLS and XMSS cannot do this. The only scalable solution is to produce a ZK Proof stating “All XMSS validator signatures are valid”. This is what LeanVM does.

Thus, many hopes are placed on ZK Proofs. However, they:

- are relatively poorly explored, still novel territory
- are very expensive to produce!

## Poseidon

The Poseidon algorithm was developed specifically with ZKP-friendliness in mind. But before it gets adopted, it needs to be rigorously tested. These initiatives test the limits of Poseidon 1 and Poseidon 2:

[Poseidon 1 Collision Prize](https://www.poseidon-initiative.info/#h.fx9618au391i)
[Bounty Program 2026: CICO Problem](https://www.poseidon-initiative.info/#h.l22wkxnie0q7)
[Bounty Program 2026: Density Problem](https://www.poseidon-initiative.info/#h.2tqzbg6s1xb5)
[Bounty Program 2026: Zero-test Problem](https://www.poseidon-initiative.info/#h.dibj06k8bvsh)

## Proximity Prize

ZK proofs rely on the Reed-Solomon Coding Scheme.

Data is represented via numbers, e.g. `[4, 9, 13]`. These numbers can be coefficients of a polynomial, e.g. $f(x)=4x^2 + 9x^1 + 13x^0$. Simply put, when a verifier verifies a ZK proof, they verify the values of the polynomial at certain points: $f(2)\neq 47 \rightarrow \text{forged proof}$.

The question is how many points we need to evaluate, in order to be certain that the proof hasn’t been forged? We can’t verify all points, because it would be too expensive!

The Reed-Solomon coding scheme evaluates a high number of points within a margin that is considered safe. But what if we can bring this safe interval way down? This would speed up zero-knowledge proofs significantly, thus making them applicable to a wide range of performance-critical problems (e.g. prove a person is an adult without exposing their date of birth).

The two problems, posed by the proximity prize, aim to provide a definitive evaluation of Reed-Solomon parameters within a safe margin.

[Proximity Prize](https://proximityprize.org/)

## Our Job

1. Study the theory behind these Poseidon and Reed-Solomon Coding Schemes
2. Try to come up with a solution to any of the open problems in the initiatives. What we choose to tackle depends on which algorithms we understand better and what ideas to “hack” them we have.
3. Even if we never reach a solution, we should document any useful findings, in order to contribute to this space.

## Resources

- [Algebra](./Algebra.md) - a short list of important algebraic concepts
