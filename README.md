# Tema 1 APD

### Paralelizare Rescale

- Se verifica in main daca imaginea trebuie scalata.
    - Daca trebuie scalata se apeleaza rescale_image in P thread-uri
    - Daca nu se trece mai departe
- In main se aloca memorie pentru noua imagine
- Inainte de crearea threadurilor se initializeaza vectorul de structuri Args ce contine argumentele functiei.
    - imaginea initiala
    - noua imagine
    - id-ul thread-ului
    - numarul de thread-uri
- In functie se paralizaza primul for cu formula: **start = ID * (double)N / P; end = min((ID + 1) * (double)N / P, N);** unde N este rezolutia pe axa OX si P numarul de thread-uri
- In main dupa ce se termina thread-urile se elibereaza zona de memorie ocupata de vechea imagine.
- Se continua algoritmul initial
