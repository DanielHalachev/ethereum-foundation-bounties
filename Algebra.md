# Algebra

A few useful concepts from algebra that we will need to understand the Poseidon algorithm and the Reed-Solomon coding scheme.

## Group

A group $G$ in algebra is a set of elements together with a binary operation $\circledast$ on G that combines two elements in $G$ to form a third element in G satisfying the following axioms:

- associativity: $(a \circledast b) \circledast c = a \circledast (b \circledast c)$
- identity element: there exists an element $e \in G$ such that $e \circledast a = a \circledast e = a \forall a \in G$
- inverse element: $\forall a \in G\space\exists\space a^{-1}\in G$ such that $a \circledast a^{-1} = a^{-1} \circledast a = e$

### Examples

- $(\mathbb{Z}, +)$, the set of integers with addition as the operation
- $(\mathbb{Z}^{*}, \cdot)$, is not a group, because there is no inverse element for integers that is also an integer (e.g. $2^{-1} = \frac{1}{2}$, which is not an integer)
- $(\mathbb{Q}^{*}, \cdot)$, the set of non-zero rational numbers with multiplication as the operation, is a group
- $(\mathbb{Q}, \cdot)$, is not a group, because there is no inverse element for zero that is also a rational number (e.g. $0^{-1} = \frac{1}{0}$, which is not a rational number)
- ($\mathbb{Z}_4, +)$, the set of integers modulo 4 with addition as the operation is a group
- $(\mathbb{Z}_4^{*}, \cdot)$, the set of integers modulo 4 with multiplication as the operation is not a group, because there is no inverse element for 2 that is also an integer modulo 4:
  - $2\circledast 1=2 \pmod{4}$
  - $2\circledast 2=0 \pmod{4}$
  - $2\circledast 3=2 \pmod{4}$
  - $2\circledast 4=0 \pmod{4}$

$2$ does not have an inverse element in $\mathbb{Z}_4^{*}$, because $2$ is not coprime to $4$.

### Groups for integers under multiplication

To form a group $(\mathbb{Z}_p^*, \cdot)$ under multiplication, $p$ needs to be a prime number.

### Generator

$g$ is a generator of a group $G$ if every element of $G$ can be obtained by repeatedly applying the group operation to $g$.

For $+$, we apply the operation $+$ $n$ times, thus we check $ng$.

For $\cdot$, we apply the operation $\cdot$ $n$ times, thus we check $g^n$.

For example, $2$ is not a generator of $(\mathbb{Z}_7^{*}, \cdot)$:

- $2^1 = 2 \pmod{7}$
- $2^2 = 4 \pmod{7}$
- $2^3 = 1 \pmod{7} \space\implies\space \text{2 has order 3 in } (\mathbb{Z}_7^{*}, \cdot)\text{, values start cycling from here}$

But $3$ is a generator of $(\mathbb{Z}_7^{*}, \cdot)$:

- $3^1 = 3 \pmod{7}$
- $3^2 = 2 \pmod{7}$
- $3^3 = 6 \pmod{7}$
- $3^4 = 4 \pmod{7}$
- $3^5 = 5 \pmod{7}$
- $3^6 = 1 \pmod{7}\space \implies \space \text{3 has order 6 in } (\mathbb{Z}_7^{*}, \cdot)\text{, values start cycling from here}$

### Abelian group

A group $G$ is called an abelian group if the group operation is commutative, i.e. $a \circledast b = b \circledast a \space\space \forall a, b \in G$. I.e., the order of the operands does not matter.

## Ring

A ring $R$ in algebra is a set of elements together with two binary operations $+$ and $\cdot$ on $R$ satisfying the following axioms:

- $(R, +)$ is an abelian group
- associativity: $\forall a, b, c \in R:\space a\cdot(b\cdot c) = (a\cdot b)\cdot c$
- distributivity: $a \cdot (b + c) = a \cdot b + a \cdot c$ and $(a + b) \cdot c = a \cdot c + b \cdot c \space\space \forall a, b, c \in R$

Thus, $R$ does not need to have inverse elements. Thus R does not need to be a group under multiplication $\cdot$, but it needs to be a group under addition $+$.

### Examples

- $(\mathbb{Z}, +, \cdot)$, the set of integers with addition and multiplication as the operations
- $(\mathbb{Z}_n, +, \cdot)$, the set of integers modulo $n$ with addition and multiplication as the operations

## Field

A field $F$ in algebra is a set of elements together with two binary operations $+$ and $\cdot$ on $F$ satisfying the following axioms:

- $(F, +)$ is an abelian group under addition
- $F^{*} = F \setminus \{0\}$ is an abelian group under multiplication.

### Examples

- $(\mathbb{Z}_p, +, \cdot)$, when $p$ is prime, is a field.

### Finite Fields

A field $F_q$ is finite, if $q$ is a prime number, or a power of a prime number, i.e. $q = p^n$, where $p$ is a prime number and $n$ is a positive integer.\
