# Proyek Simulasi H2: Variational Monte Carlo (VMC) & Variational Quantum Eigensolver (VQE)

Repositori ini berisi implementasi dan perbandingan metode klasik dan kuantum untuk menghitung energi ground state dari molekul Hidrogen ($H_2$).

## Struktur File dan Penjelasan

Proyek ini terdiri dari dua Jupyter Notebook utama:

### 1. `H2Classic.ipynb` (Klasik VMC)
Notebook ini mengimplementasikan metode **Variational Monte Carlo (VMC)** secara klasik untuk molekul $H_2$.

*   **Fokus Utama**: Menghitung energi menggunakan fungsi gelombang trial dengan parameter variasional $\alpha$.
*   **Metode Optimasi**:
    *   **Grid Search**: Memindai nilai $\alpha$ secara manual untuk menemukan energi minimum.
    *   **Gradient Descent (On-the-fly)**: Mengupdate $\alpha$ secara iteratif berdasarkan gradien energi.
*   **Fitur**:
    *   Perhitungan energi lokal dan fungsi gelombang.
    *   Analisis grafik "Energy vs Trial Step" untuk melihat stabilitas simulasi.
    *   Perbandingan hasil VMC (~ -1.127 Hartree) dengan nilai eksak (~ -1.174 Hartree), di mana selisih disebabkan oleh kurangnya korelasi elektron-elektron pada ansatz yang digunakan.

### 2. `H2VQE.ipynb` (Kuantum VQE)
Notebook ini mengimplementasikan **Variational Quantum Eigensolver (VQE)** menggunakan framework **Qiskit Primitives V2**.

*   **Sistem**: Hamiltonian 1-qubit untuk $H_2$ pada jarak ikatan 0.735 Å.
*   **Ansatz**: Sirkuit parameterisasi 3-parameter: $R_x(\theta_0) R_z(\theta_1) R_x(\theta_2)$.
*   **Estimator**: Menggunakan `StatevectorEstimator` untuk simulasi eksak (tanpa shot noise).
*   **Fitur**:
    *   **Custom VQE Logger**: Callback function untuk memantau iterasi optimasi, nilai energi, dan perubahan parameter secara real-time.
    *   **Visualisasi**: Grafik konvergensi optimasi dan analisis landscape energi.
*   **Hasil**: Mampu mencapai energi ground state yang sangat akurat (-1.145977 Hartree untuk Hamiltonian yang direduksi).

## Dependensi

Proyek ini membutuhkan pustaka Python berikut:
*   `numpy`
*   `scipy`
*   `matplotlib`
*   `qiskit`

## Cara Penggunaan

1.  Buka notebook yang diinginkan (`H2Classic.ipynb` untuk simulasi klasik atau `H2VQE.ipynb` untuk simulasi kuantum).
2.  Jalankan sel secara berurutan.
3.  Untuk VQE, perhatikan output log yang menunjukkan penurunan energi di setiap iterasi optimizer.
