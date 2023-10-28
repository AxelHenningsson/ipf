import setuptools

with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setuptools.setup(
    name="ipf",
    version="0.0.1",
    author="Axel Henningsson",
    author_email="nilsaxelhenningsson@gmail.com",
    description="A bare-bones implementation of the inverse pole figure.",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/AxelHenningsson/ipf",
    project_urls={
        "Documentation": "https://github.com/AxelHenningsson/ipf",
    },
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.8,<3.11",
    install_requires=["numpy",
                      "matplotlib",
                      "shapely",
                      "scipy",
                      "xfab"
                      ]
)
