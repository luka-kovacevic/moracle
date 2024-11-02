from setuptools import setup, find_packages

setup(name="moracle",  # Replace with your package name
    version="0.1.0",
    description="Molecule Oracle (MOracle): A tool for molecule-protein binding affinity analysis",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/yourusername/moracle",  # Replace with the URL of your repository
    packages=find_packages("moracle"),  # Assumes your code is in src/; adjust if necessary
    package_dir={"": "moracle"},
    include_package_data=True,
    install_requires=[
        "streamlit>=1.0",
        "pandas>=1.0",
        "numpy>=1.18.0",
        "streamlit_ketcher",  # Include other dependencies as needed
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)