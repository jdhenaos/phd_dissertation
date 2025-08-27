build-invitro:
	cd invitro && Rscript -e "workflowr::wflow_build()"
	mkdir -p docs/invitro
	rsync -av invitro/docs/ docs/invitro/

