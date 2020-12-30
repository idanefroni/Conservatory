%define _prefix /usr

URL: http://bioweb.me/maffilter

Name: maffilter
Version: 1.3.1
Release: 1%{?dist}
License: CECILL-2.0
Vendor: The Bio++ Project
Source: %{name}-%{version}.tar.gz
Summary: The Multiple Alignment Format file processor
Group: Productivity/Scientific/Other

Requires: libbpp-phyl-omics3 = 2.4.1
Requires: libbpp-seq-omics3 = 2.4.1
Requires: libbpp-phyl12 = 2.4.1
Requires: libbpp-seq12 = 2.4.1
Requires: libbpp-core4 = 2.4.1
Requires: zlib

BuildRoot: %{_builddir}/%{name}-root
BuildRequires: cmake >= 2.8.11
BuildRequires: gcc-c++ >= 4.7.0
BuildRequires: groff
BuildRequires: texinfo >= 4.0.0
BuildRequires: libbpp-core4 = 2.4.1
BuildRequires: libbpp-core-devel = 2.4.1
BuildRequires: libbpp-seq12 = 2.4.1
BuildRequires: libbpp-seq-devel = 2.4.1
BuildRequires: libbpp-phyl12 = 2.4.1
BuildRequires: libbpp-phyl-devel = 2.4.1
BuildRequires: libbpp-seq-omics3 = 2.4.1
BuildRequires: libbpp-seq-omics-devel = 2.4.1
BuildRequires: libbpp-phyl-omics3 = 2.4.1
BuildRequires: libbpp-phyl-omics-devel = 2.4.1
BuildRequires: zlib-devel

%if 0%{?fedora} >= 22 || 0%{?mageia} <= 6
Requires: bzip2-libs
BuildRequires: bzip2-devel
%else
Requires: libbz2
BuildRequires: libbz2-devel
%endif

AutoReq: yes
AutoProv: yes
%if 0%{?mdkversion}
%if 0%{?mdkversion} >= 201100
BuildRequires: xz
%define zipext xz
%else
BuildRequires: lzma
%define zipext lzma
%endif
%else
BuildRequires: gzip
%define zipext gz
%endif

#Distribution specific requirements:

%if 0%{?suse_version} >= 1200 && 0%{?suse_version} < 1300 
Requires: libboost_iostreams1_49_0 
BuildRequires: boost-devel
%else
%if 0%{?suse_version} >= 1300
Requires: libboost_iostreams1_53_0 
BuildRequires: boost-devel
%else
Requires: libboost-iostreams 
BuildRequires: boost-devel
%endif
%endif


%description
MafFilter is a program for processing and analysing files in the Multiple Alignment Format (MAF).
A description of MAF files can be found on the UCSC genome browser (http://genome.ucsc.edu/FAQ/FAQformat.html#format5).
MafFilter can be used to design a pipeline as a series of consecutive filters, each performing a dedicated analysis.
Many filters are available, from alignment cleaning to phylogeny reconstruction and population genetics analysis.

%prep
%setup -q

%build
CFLAGS="$RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix} -DCOMPRESS_PROGRAM=%{compress_program}"
cmake $CMAKE_FLAGS .
make

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root)
%doc AUTHORS LICENSE INSTALL ChangeLog
%{_prefix}/bin/maffilter
%{_prefix}/share/man/man1/maffilter.1.%{zipext}
%{_prefix}/share/info/maffilter.info.%{zipext}

%changelog
* Fri Aug 17 2018 Julien Dutheil <dutheil@evolbio.mpg.de> 1.3.1-1
- Compatibility update with bio++ 2.4.1
- Documentation update
* Wed Mar 14 2018 Julien Dutheil <dutheil@evolbio.mpg.de> 1.3.0-1
* Fri Jun 09 2017 Julien Dutheil <dutheil@evolbio.mpg.de> 1.2.1-1
* Wed May 24 2017 Julien Dutheil <dutheil@evolbio.mpg.de> 1.2.0-1
* Fri Sep 26 2014 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.1.0-1
- Initial spec file.
